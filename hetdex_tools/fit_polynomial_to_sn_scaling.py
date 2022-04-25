"""

Fit a polynomial of S/N cut to
noise multipliers

"""
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 10.0
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
from numpy import polyfit, polyval, array



def fit_karls_old_scaling_v1():
    """
    Fit completeness files like this

    #name  s/n ccut0 ccut1  f50x
    sn4.8.use 4.8 0.000 0.0000 4.49 4.607
    sn5.0.use 5.0 0.000 0.0000 4.60 4.783
    sn5.5.use 5.5 0.000 0.0000 4.82 5.282
    sn6.0.use 6.0 0.000 0.0000 5.18 5.805
    sn6.5.use 6.5 0.000 0.0000 5.61 6.362
    sn7.0.use 7.0 0.000 0.0000 6.04 6.963

    """

    sn = []
    scaling = []
    with open("karl_file", 'r') as fp:
        for line in fp:
            els = line.strip().split()
            if "#" not in els[0]:
                sn.append(float(els[1]))
                scaling.append(float(els[4]))
    
    sn = array(sn)
    scaling = array(scaling)
    
    p = polyfit(sn, scaling, 4)
    print(p)
    for s, scl in zip(sn, scaling):
        print("{:2.1f} {:4.3f} {:4.3f} {:4.3f}".format(s, scl, polyval(p, s), scl/polyval(p, s)))
    
    
    plt.plot(sn, scaling, "k*")
    plt.plot(sn, polyval(p, sn), "r--")
    plt.ylabel("Noise scaling")
    plt.xlabel("S/N cut")
    plt.tight_layout()
    plt.show()

def fit_karl_wavelength_dependence():
    """
    Fit 1/f_w, the wavelength dependence
    derived from Karl's simulations


    *old version (v2)*

    3500 1.105
    3700 1.073
    3950 1.046
    4200 1.032
    4500 1.032
    4750 1.023
    5050 1.020
    5300 1.022
    5500 1.020

    *new version (v4)*
    April 13 2022
    Note: definition changed was
    1/ this value before

    3500 0.98
    3957 0.85
    4229 0.81
    4500 0.79
    4771 0.79
    5043 0.80
    5500 0.83

    """

    lambda_ = array([3500, 3957, 4229, 4500, 4771, 5043, 5500])
    f_karl = array([0.98, 0.85, 0.81, 0.79, 0.79, 0.80, 0.83])

    p = polyfit(lambda_, f_karl, 3)
    print(f"Resulting polynomial terms {p}")

    for w, f in zip(lambda_, f_karl):
        print("{:2.1f} {:4.3f} {:4.3f} {:4.3f}".format(w, f, polyval(p, w), f/polyval(p, w)))
    
    # f_karl is f_lambda
    scale = 1.5
    plt.figure(figsize=(3.5*scale,3*scale))
    plt.plot(lambda_, f_karl, "k*", 
             label="Source Simulations")
    plt.plot(lambda_, polyval(p, lambda_), "r--",
             label="Polynomial fit")
    plt.ylabel("$F_{\lambda}(\lambda)$", fontsize=10.0)
    plt.xlabel("Wavelength $\lambda$ (Angstrom)", fontsize=10.0)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    fit_karl_wavelength_dependence()
    #fit_karls_old_scaling_v1()
