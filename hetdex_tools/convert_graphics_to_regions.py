import numpy as np
from regions import PixCoord, PolygonSkyRegion, PolygonPixelRegion, RectangleSkyRegion
from astropy.coordinates import SkyCoord

def read_graphics_file(infile):
    # read all the points from the graphics file

    f = open(infile, "r")
    
    for i, line in enumerate(f):
        if i == 0:
            numpoly = line.split()[0]
            print("numpoly = {}".format(numpoly))
            # parse polygons
            polygons_ra = []
            polygons_dec = []
            continue
            
        if i <= 1:
            continue
            
        s = np.array( line.split())
        #print(s)
        # skip "graphics" lines, for now
        if s[0] == "graphics":
            numpoints = int(s[3])
            #print(numpoints)
            continue

        # save the points
        ra = []
        dec = []
        for m in np.arange(0, numpoints):
            ra.append( float(s[2*m])) 
            dec.append( float(s[2*m+1]))

        polygons_ra.append(ra)
        polygons_dec.append(dec)

    return polygons_ra, polygons_dec

graphics_file = '/work2/00436/djeong/stampede2/decodeHETDEX/HETDEX/DR4/IFUdata/polygons/cosmos_hdr4_mask-raw.graphics'

polygons_ra, polygons_dec = read_graphics_file(graphics_file)
hetdex_regions = []

for poly_ra, poly_dec in zip(polygons_ra, polygons_dec):
    hetdex_regions.append(PolygonSkyRegion(vertices=SkyCoord(poly_ra, poly_dec, unit='deg', frame='fk5')))
