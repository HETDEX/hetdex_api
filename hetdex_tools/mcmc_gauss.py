"""
Cloned from ELiXer mcmc_gauss.py

Fits an emission line with a simple, single Gaussian for four parameters:
line center (mu), line width (sigma), area (A), and y offset (y)
    * Note the y offset is not a line fit (that is, the slope is zero)

"""

from __future__ import print_function

import numpy as np
from scipy.signal import medfilt
import emcee
import copy
import warnings
import logging
from hetdex_api.input_utils import setup_logging

#uncomment the next three includes AND uncomment the visualize function at the bottom if you want them
#import matplotlib.pyplot as plt
#import corner
#import io

SNR_SIGMA_WIDTH = 2 #line center +/- this width over which to compute the noise
SNR_LINEFLUX_FRACTION = 0.954 #0.682 for 1 sigma, 0.954 for 2 sigma, 0.996 for 3 sigma, assume 1.0 thereafter
UncertaintyRange = [16,50,84] #e.g. for the uncertainty as distribution

def getnearpos(array,value):
    """
    Nearest, but works best (with less than and greater than) if monotonically increasing. Otherwise,
    lt and gt are (almost) meaningless

    :param array:
    :param value:
    :return: nearest index, nearest index less than the value, nearest index greater than the value
            None if there is no less than or greater than
    """
    if type(array) == list:
        array = np.array(array)

    idx = (np.abs(array-value)).argmin()

    if array[idx] == value:
        lt = idx
        gt = idx
    elif array[idx] < value:
        lt = idx
        gt = idx + 1
    else:
        lt = idx - 1
        gt = idx

    if lt < 0:
        lt = None

    if gt > len(array) -1:
        gt = None

    return idx, lt, gt

class MCMC_Gauss:

    def __init__(self,logger=None):
        #intial values are meant to be near the truth
        #and are expected to come from, say, some initial "best" fit


        self.log = logger
        if self.log is None:
            self.log = setup_logging()
            self.log.setLevel(logging.INFO)

        self.initial_mu = None
        self.initial_sigma = None
        self.initial_A = None  #note: set to a negative value if this is an absorption line
        self.initial_y = None
        self.initial_peak = None

        self.max_sigma = 20.0
        self.range_mu = 5.0
        self.max_A_mult = 2.0
        self.max_y_mult = 2.0
        self.min_y = -10.0

        self.data_x = None
        self.data_y = None
        self.err_x = None
        self.err_y = None

        #just for reference ... MCMC itself does not need to know about this
        #the caller DOES though and needs to adjust the line_flux accordingly
        #self.dx = None #original bin width IF NOT part of the data_y already

        #this is mostly a guess ... no great way to automate, but this is pretty quick
        #and since the initials are from a scipy curve fit, we stabablize pretty fast
        self.burn_in = 100
        self.main_run = 1000
        self.walkers = 100

        self.sampler = None #mcmc sampler
        self.samples = None #resulting samples

        #####################
        # Outputs
        #####################
        #3-tuples [0] = fit, [1] = fit +16%,  [2] = fit - 16%
        self.mcmc_mu = None
        self.mcmc_sigma = None
        self.mcmc_A = None  #note: typically this is over the HETDEX 2AA bins, so if using as the area as integrated
                            # lineflux you need to divide by 2AA and scale appropriately (e.g. as 1e-17)
        self.mcmc_y = None

        #not tuple, just single floats
        self.mcmc_snr = None
        self.mcmc_snr_err = 0


    def approx_symmetric_error(self,parm): #parm is assumed to be a 3 vector as [0] = mean, [1] = +error, [2] = -error

        try:
            if parm is None or (len(parm)!= 3) :
                return None

            p1 = abs(parm[1])
            p2 = abs(parm[2])
            avg = 0.5*(p1+p2)

            if avg == 0:
                return 0

            similarity = abs(p1-p2)/avg

            if similarity > 0.1:
                self.log.warning("Warning! Asymmetric uncertainty similarity = %0.3g (%0.3g, %0.3g)" %(similarity,p1,p2))

            #for now, do it anyway
            return avg
        except:
            return None


    def noise_model(self):
        #todo: fill in some specialized model for the noise
        return 0.0

    def compute_model(self,x,mu, sigma, A, y):
        try:
            return A * (np.exp(-np.power((x - mu) / sigma, 2.) / 2.) / np.sqrt(2 * np.pi * sigma ** 2)) + y
        except:
            return np.nan

    def model(self,x,theta):
        mu, sigma, A, y, ln_f = theta #note: not using ln_f here
        if (x is None) or (mu is None) or (sigma is None):
            return None
        try:
            value = self.compute_model(x,mu, sigma, A, y)
            # note: noise is separate and included in the lnlike() function
        except:
            value = np.nan
        return value

    def lnlike(self, theta, x, y, yerr):
        ln_f = theta[-1] #last parameter in theta
        model = self.model(x, theta)
        noise = self.noise_model()
        diff = y - (model + noise)
        #assumes some additional uncertainties in y based on an underestimation in the model by some factor f
        # assume that the (distribution of) errors in y are known and independent
        sigma2 = (self.err_y ** 2)
        return -0.5 * (np.sum((diff ** 2) / sigma2 + np.log(sigma2)))

    # if any are zero, the whole prior is zero
    # all priors here are uniformitive ... i.e they are all flat ... either zero or one
    def lnprior(self, theta):  # theta is a n-tuple (_,_,_ ... )
        mu, sigma, A, y, ln_f = theta
        # note: could take some other dynamic maximum for y (like compute the peak ... y can't be greater than that

        if self.initial_A < 0 : #same as emission, but "A" is negative (flip sign) and y is between a max and zero
            if (-self.range_mu < mu - self.initial_mu < self.range_mu) and \
                    (0.0 < sigma < self.max_sigma) and \
                    (self.max_A_mult * self.initial_A < A < 0.0) and \
                    (self.min_y < y < self.max_y_mult * self.initial_peak):
                return 0.0  # remember this is ln(prior) so a return of 0.0 == 1  (since ln(1) == 0.0)
        else:
            if (-self.range_mu < mu - self.initial_mu < self.range_mu) and \
                    (0.0 < sigma < self.max_sigma) and \
                    (0.0 < A < self.max_A_mult * self.initial_A) and \
                    (self.min_y < y < self.max_y_mult * self.initial_peak):
                return 0.0  # remember this is ln(prior) so a return of 0.0 == 1  (since ln(1) == 0.0)
        return -np.inf  # -999999999 #-np.inf #roughly ln(0) == -inf

    def lnprob(self, theta, x, y, yerr):
        """
        ln(probability)

        :param theta: parameters to check
        :param x:  THE data (x axis or wavelengths, in this case)
        :param y: THE data (y axis or flux counts, in this case)
        :param yerr:  The error on the y axis data flux counts
        :return:
        """
        lp = self.lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.lnlike(theta, x, y, yerr)  # again, since logs, this is a sum ln(likelihood) + ln(prior)


    def sanity_check_init(self):
        """
        evaluation that the initialization data makes sense before running MCMC
        :return:
        """
        try:
            #if the error on y is None or if it is all zeros, set to all ones
            if self.err_y is None:
                self.err_y = np.ones(np.shape(self.data_y))
            elif not np.any(self.err_y):
                self.err_y = np.ones(np.shape(self.data_y))

            if self.err_x is None:
                self.err_x = np.ones(np.shape(self.data_x))

            if (self.data_x is None) or (self.data_y is None) or (self.err_y is None):
                self.log.debug("Sanity check failed. data_x or data_y or err_y is None")
                return False

            if len(self.data_x) == len(self.data_y) == len(self.err_y): #leave off self.err_x as could be None
                if (self.err_x is not None):
                    if len(self.data_x) != len(self.err_x):
                        self.log.debug("Sanity check failed. len(data_x) != len(err_x)")
                        return False

                if (self.initial_sigma is None) or (self.initial_mu is None) or (self.initial_A is None) or (self.initial_y is None):
                    self.log.debug("Sanity check failed. initial sigma, mu, A, or y is None")
                    return False

                if self.initial_peak is None:
                    left,*_ = getnearpos(self.data_x,self.initial_mu-self.initial_sigma*4)
                    right,*_ = getnearpos(self.data_x,self.initial_mu+self.initial_sigma*4)

                    self.initial_peak = max(self.data_y[left:right])

                if self.initial_sigma < 0.0: #self.initial_A < 0.0  ... actually, leave A alone .. might allow absorportion later
                    self.log.debug("Sanity check failed. initial sigma < 0")
                    return False

                if ((self.initial_A > 0) and (self.initial_y > self.initial_peak) or \
                        (self.initial_A < 0) and (self.initial_y < self.initial_peak) ):
                    #i.e. if an emission (A > 0) then y must be less than the peak
                    # else if an absorption line (A < 0) then y must be greater than the peak
                    self.log.debug("Sanity check failed. y offset inconsistent with area (A) and initial peak.")
                    return False
            else:
                self.log.debug("Sanity check failed. lengths of data_x, data_y, and err_y do not match.")
                return False
            return True
        except:
            self.log.warning("Exception in mcmc_gauss::sanity_check",exc_info=True)
            return False

    def run_mcmc(self):

        #cannot have nans
        #note: assumes data_x (the spectral axis) and err_x have none since they are on a known grid
        data_nans = np.isnan(self.data_y)
        err_nans = np.isnan(self.err_y)

        if (np.sum(data_nans) > 0) or (np.sum(err_nans) > 0):
            self.data_y = copy.copy(self.data_y)[~data_nans]
            self.err_y = copy.copy(self.err_y)[~data_nans]
            self.data_x = copy.copy(self.data_x)[~data_nans]
            self.err_x = copy.copy(self.err_x)[~data_nans]
            #and clean up any other nan's in the error array for y
            err_nans = np.isnan(self.err_y)
            self.err_y[err_nans] = np.nanmax(self.err_y*10)

        if not self.sanity_check_init():
            self.log.info("Sanity check failed. Cannot conduct MCMC.")
            return False

        result = True

        #here for initial positions of the walkers, sample from narrow gaussian (hence the randn or randNormal)
        #centered on each of the maximum likelihood selected parameter values
        #mu, sigma, A, y, ln_f = theta #note f or ln(f) is another uncertainty ...an underestimation of the variance
        #                               by some factor (f) .... e.g. variance = variance + f * model
        initial_pos = [self.initial_mu,self.initial_sigma,self.initial_A,self.initial_y,0.0]

        #even with the random nudging the pos values must be greater than (or less than for absorption) these values

        #mostly for the A (area)
        if self.initial_A < 0: #absorber
            max_pos = [np.inf, np.inf,     0.0, max(self.data_y),  np.inf]
            min_pos = [   0.0,   0.01, -np.inf,          -np.inf, -np.inf]
        else:
            #here, because of the max check, none mu, sigma, or A will be negative
            max_pos = [np.inf, np.inf,np.inf,max(self.data_y), np.inf] #must be less than this
            min_pos = [   0.0,  0.01,   0.01,         -np.inf,-np.inf] #must be greater than this

        ndim = len(initial_pos)
        scale = np.array([10.,5.,2.0*self.initial_A,5.0*self.initial_y,0.01]) #don't nudge ln_f ...note ln_f = -4.5 --> f ~ 0.01

        try:
            ##############################################################################
            # This is an alternate way to control the jitter in the initial positions,
            # can uncomment this block if the jitter is not sufficient
            ##############################################################################
            # ip_mu = initial_pos[0] + np.random.uniform(-1.0*(self.data_x[1]-self.data_x[0]),1.0*(self.data_x[1]-self.data_x[0]),self.walkers)
            # #sigma cannot go zero or below
            # ip_sigma = initial_pos[1] + np.random.uniform(-0.5*self.initial_sigma,0.5*self.initial_sigma,self.walkers)
            # #area cannot flip signs
            # ip_A = initial_pos[2] +  np.random.uniform(0,1.0*self.initial_A,self.walkers)
            # #y should not exceed min/max data value, but won't cause and error if it does
            # # ... should technically never be negative regardless of absorption or emission
            # ip_y = np.random.uniform(0,max(self.data_y),self.walkers)
            # ip_lnf = np.zeros(self.walkers) #np.random.uniform(0.005,0.015,self.walkers) #np.zeros(self.walkers)
            #
            # #for p in pos: #just a debug check
            # #  print(f"{p[0]:0.4g},{p[1]:0.4g},{p[2]:0.4g},{p[3]:0.4g},{p[4]:0.4g}")

            ##############################################################################
            # OTHERWISE, keep the line below
            #
            ##############################################################################

            pos = [np.minimum(np.maximum(initial_pos + scale * np.random.uniform(-1,1,ndim),min_pos),max_pos) for i in range(self.walkers)]

            #build the sampler
            self.sampler = emcee.EnsembleSampler(self.walkers, ndim, self.lnprob,
                                                 args=(self.data_x,self.data_y, self.err_y))
            #args are the positional args AFTER theta for self.lnprob function

            with warnings.catch_warnings(): #ignore the occassional warnings from the walkers (NaNs, etc that reject step)
                warnings.simplefilter("ignore")
                self.log.debug("MCMC burn in (%d) ...." %self.burn_in)
                pos, prob, state = self.sampler.run_mcmc(pos, self.burn_in,skip_initial_state_check=False)  # burn in
                self.log.debug("MCMC main run (%d) ..." %self.main_run)
                pos, prob, state = self.sampler.run_mcmc(pos, self.main_run, rstate0=state,skip_initial_state_check=False)  # start from end position of burn-in

            self.samples = self.sampler.flatchain  # collapse the walkers and interations (aka steps or epochs)

            self.log.debug("MCMC mean acceptance fraction: %0.3f" %(np.mean(self.sampler.acceptance_fraction)))

            #for each, in order
            #v[0] is the 16 percentile (~ - 1sigma)
            #v[1] is the 50 percentile (so the "average")
            #v[2] is the 84 percentile (~ +1sigma)
            #the tuple reports then as ["average", "84th - average", "average - 16th"]
            #should always be positive (assuming a positive value) BUT when printed to the log, the 3rd value is made
            #to be negative showing that you would subtract it from the average to get to the 16th percentile

            #using 68% interval
            self.mcmc_mu, self.mcmc_sigma, self.mcmc_A, self.mcmc_y, mcmc_f = \
                map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),zip(*np.percentile(self.samples, UncertaintyRange,axis=0)))

            try: #basic info used by multiple SNR calculations
                bin_width = self.data_x[1] - self.data_x[0]
                left,*_ = getnearpos(self.data_x,self.mcmc_mu[0]-self.mcmc_sigma[0]*SNR_SIGMA_WIDTH)
                right,*_ = getnearpos(self.data_x,self.mcmc_mu[0]+self.mcmc_sigma[0]*SNR_SIGMA_WIDTH)

                #at 4 sigma the mcmc_A[0] is almost identical to the model_fit (as you would expect)
                #note: if choose to sum over model fit, remember that this is usually over 2AA wide bins, so to
                #compare to the error data, need to multiply the model_sum by the bin width (2AA)
                #(or noting that the Area == integrated flux x binwidth)

                #model_fit = self.compute_model(self.data_x[left:right],self.mcmc_mu[0],self.mcmc_sigma[0],self.mcmc_A[0],self.mcmc_y[0])
                data_err = copy.copy(self.err_y[left:right])
                data_err[data_err<=0] = np.nan #Karl has 0 value meaning it is flagged and should be skipped

                self.mcmc_snr = SNR_LINEFLUX_FRACTION*abs(self.mcmc_A[0]/bin_width) / np.sqrt(np.nansum(data_err**2))
                self.mcmc_snr_err = abs(0.5*(self.mcmc_A[1]+self.mcmc_A[2])/self.mcmc_A[0] * self.mcmc_snr)
                self.log.info(f"MCMC SNR model Area with data error: {self.mcmc_snr} +/- {self.mcmc_snr_err}")

            except:
                self.log.warning("Exception calculating MCMC SNR: ", exc_info=True)

            if self.mcmc_snr is None:
                self.mcmc_snr = -1

            #note: these are logged as ["avg", +err, -err] so the last value becomes the negative
            self.log.info("MCMC mu: initial[%0.5g] mcmc(%0.5g, +%0.5g, -%0.5g)" %
                     (self.initial_mu, self.mcmc_mu[0],self.mcmc_mu[1],self.mcmc_mu[2]))
            self.log.info("MCMC sigma: initial[%0.5g] mcmc(%0.5g, +%0.5g, -%0.5g)" %
                     (self.initial_sigma, self.mcmc_sigma[0],self.mcmc_sigma[1],self.mcmc_sigma[2]))
            self.log.info("MCMC A: initial[%0.5g] mcmc(%0.5g, +%0.5g, -%0.5g) *usually over 2AA bins" %
                     (self.initial_A, self.mcmc_A[0],self.mcmc_A[1],self.mcmc_A[2] ))
            self.log.info("MCMC y: initial[%0.5g] mcmc(%0.5g, +%0.5g, -%0.5g)"%
                     (self.initial_y, self.mcmc_y[0],self.mcmc_y[1],self.mcmc_y[2]))
            self.log.info("MCMC SNR: %0.5g" % self.mcmc_snr)
            self.log.info("MCMC f: initial[%0.5g] mcmc(%0.5g, +%0.5g, -%0.5g)" %
                     (0.0, mcmc_f[0], mcmc_f[1], mcmc_f[2]))
        except:
            self.log.error("Exception in mcmc_gauss::run_mcmc",exc_info=True)
            result = False

        return result

    #need to uncomment matplotlib, corner, and io at the top if you want this function
    # def visualize(self,filename=None):
    #     try:
    #         if self.samples is not None:
    #             warnings.simplefilter(action='ignore', category=FutureWarning)
    #
    #             fig = corner.corner(self.samples, labels=["$mu$", "$sigma$", "$A$", "$y$","f"],
    #                                 truths=[self.initial_mu, self.initial_sigma, self.initial_A, self.initial_y,None])
    #             #fifth = None is for the 'f' parameter ... there is no initial for it
    #             if filename is not None:
    #                 self.log.info('Writing: ' + filename)
    #                 fig.savefig(filename)
    #             else:
    #                 plt.show()
    #
    #             buf = None
    #             try:
    #                 buf = io.BytesIO()
    #                 fig.savefig(buf, format='png', dpi=300)
    #             except:
    #                 self.log.warning("Exception in mcmc_gauss::visualize",exc_info=True)
    #             return buf
    #     except:
    #         self.log.warning("Exception in mcmc_gauss::visualize",exc_info=True)
    #         return None
