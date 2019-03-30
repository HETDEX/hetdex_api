# -*- coding: utf-8 -*-
"""

Script to read in exlier catalogs

Created: 2019/03/28


"""

class Detection:
    #represents a single detection (may or may not be an acceptable LAE
    def __init__(self):

        self.detectid = None
        self.entryid = None #unique for a run?
        self.pdfname = None
        self.detectname = None #i.e. 20180123v009_137
        self.ra = None
        self.dec = None
        self.w = None #observed wavelength

        self.z = None #assumes LyA

        self.ew_obs = None
        self.ew_rest = None #assuming LyA z
        self.line_flux = None

        self.sigma = None
        self.chi2 = None
        self.continuum = None

        self.plae_poii_hetdex = None #hetdex only data
        self.plae_poii_aperture = None
        self.aperture_mag = None
        self.aperture_filter = None

        self.neighbors = [] #empty list (of Neighbors

    @property
    def num_neighbors(self):
        try:
            return len(self.neighbors)
        except:
            return -1


class Neighbor:
    #a catalog neighbor
    def __init__(self):
        self.ra = None
        self.dec = None
        self.mag = None
        self.filter = None #filter for the magnitude
        self.dist = None #distance in arcsecs
        self.plae_poii_cat = None #P(LAE)/P(OII)


def read_elixer_catalogs(fibfn, catfn):
    """
    read the _fib.txt and _cat.txt catalogs
    create Detections for each (combining across and within catalogs as needed)
    :return: list of Detections
    Author: Dustin Davis
    """

    detections = []
    fib_version = None #can switch on version number if necessary
    cat_version = None

    #should be one line per detection in fib.txt
    #at least one line but may be more in cat.txt

    #read fib.txt first and build out list of detections
    with open(fibfn,"r") as f:
        for line in f:
            if line[0] == '#': #except get # version 1.5.0a16
                if (fib_version is None) and ('version' in line): #1st line should be the version
                    toks = line.split()
                    if (toks is not None) and (len(toks)== 3):
                        fib_version = toks[2]
                continue

            if len(line) < 100:
                continue

            toks = line.split() #white space split

            if len(toks) < 29: #must be at least 29 for _fib.txt
                continue

            d = Detection()

            d.pdfname = toks[0]
            d.detectid = np.int64(toks[1])
            d.entryid = toks[2]
            d.ra = float(toks[3])
            d.dec = float(toks[4])
            d.w = float(toks[5])
            d.sigma = float(toks[6])
            d.chi2 = float(toks[7])
            d.line_flux = float(toks[10])
            d.continuum = float(toks[12])

            d.plae_poii_hetdex = float(toks[14])

            detections.append(d)


    #then read cat.txt and match up to existing entry in detections list
    with open(catfn,"r") as f:
        for line in f:
            if line[0] == '#':
                if (cat_version is None) and ('version' in line): #1st line should be the version
                    toks = line.split()
                    if (toks is not None) and (len(toks)== 3):
                        cat_version = toks[2]
                continue

            if len(line) < 100:
                continue

            toks = line.split() #white space split

            if len(toks) < 29: #must be at least 29 for _fib.txt
                continue

            
            pdfname = toks[0]
            detectid = np.int64(toks[1])
                
            try:
                entryid = toks[2]
                
                ra = float(toks[3])
                dec = float(toks[4])
                
                num_cat_matches = int(toks[5])
                
                w = float(toks[6])
                
                match_ra = float(toks[24])
                match_dec = float(toks[25])
                dist = float(toks[26])
                mag = toks[28]
                filter = toks[29]
                if (toks[32] is not None) and (toks[32].lower() != 'none'): #can be None if could not get an aperture
                    plae_poii = float(toks[32])
                else:
                    plae_poii = -1 
                    
                #find your mates:
                #there may be several ... the 1st is the aperture
                #each additional one is a catalog match

                for m in detections:
                    if (m.entryid == entryid) and (m.detectid == detectid) and (m.ra == ra) and (m.dec == dec) and (m.w == w):
                        #match
                        if (match_ra == 666) and (match_dec == 666) and (dist == 0.0): #this is the aperture entry
                            m.plae_poii_aperture = plae_poii
                            m.aperture_filter = filter
                            m.aperture_mag = mag
                        else: #this is a matched catalog object
                            n = Neighbor()
                            n.ra = match_ra
                            n.dec = match_dec
                            n.mag = mag
                            n.filter = filter  # filter for the magnitude
                            n.dist = dist  # distance in arcsecs
                            n.plae_poii = plae_poii  # P(LAE)/P(OII)
                            
                            m.num_neighbors += 1
                            m.neighbors.append(n)
            except:
                print('Could not ingest cat info for %s' % detectid) 

    #Now we have consumed the catalog files
    return detections
