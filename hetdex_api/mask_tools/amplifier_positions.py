"""

Set or x,y rectangle corners (clockwise ordered) that
define the amplifier positions in IFU coordinates. Used
Greg's nice diagram to help derive these shapes from 

https://github.com/grzeimann/Panacea/blob/master/images/spectrograph_layout_useful_visual.pdf

I stretch out the edges of the amplifiers to +/- 30 arcseconds in
order to account for sources on the edge


.. moduleauthor:: Daniel Farrow <dfarrow@mpe.mpg.de>

"""


# List of swapped amps, indexed by IFUSLOT + AMP and
# the amplifier that would be in their position in 
# a normal IFU
swapped_around_amps = { 
                        "095RU" : "LL",
                        "095RL" : "LU",
                        "095LL" : "RL",
                        "095LU" : "RU",
                        "046RU" : "RL",
                        "046RL" : "RU" 
                      }             

# Amplifier positions. I always extend
# out to 30 arcseconds as sometimes you
# can get detections just off the IFU edge
#
amp_corners = {
               "RU" : [
                       # 1st rectangle
                       [
                       # x coords
                        [-30.0, -30.0, -12.71, -12.71],
                       # y coords
                        [14.32, 30.0, 30.0, 14.32]
                       ],
                       # 2nd rectangle ...
                       [
                        [-12.71, -12.71, 30.0, 30.0],
                        [12.12, 30.0, 30.0, 12.12]
                       ]
                      ],

               "RL" : [ 
                       [
                        [-30.0, -30.0, -12.71, -12.71],
                        [1.1, 14.32, 14.32, 1.1],
                       ],
                       [
                        [-12.71, -12.71, 0.0, 0.0],
                        [1.1, 12.12, 12.12, 1.1]
                       ],
                       [
                        [0, 0, 30.0, 30.0],
                        [-1.1, 12.12, 12.12, -1.1]
                       ]
                      ],

               "LL" : [ 
                       [
                        [-30.0, -30.0, 0, 0],
                        [-12.12, 1.1, 1.1, -12.12]
                       ],
                       [
                        [0, 0, 12.71, 12.71],
                        [-12.12, -1.1, -1.1, -12.12],
                       ],
                       [
                        [12.71, 12.71, 30.0, 30.0],
                        [-14.32, -1.1, -1.1, -14.32]
                       ]
                      ],
              
               "LU" : [
                       [
                        [-30.0, -30.0, 12.71, 12.71],
                        [-30.0, -12.12, -12.12, -30.0]
                       ],
                       [ 
                        [12.71, 12.71, 30.0, 30.0],
                        [-30.0, -14.32, -14.32, -30.0]
                       ]
                      ]
              }


amp_fiber_ranges = {
                    "RU" : [337, 448],
                    "RL" : [225, 336],
                    "LL" : [113, 224],
                    "LU" : [1, 112]
                   }

if __name__ == "__main__":

    import six
    from numpy import array
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from  pyhetdex.het.ifu_centers import IFUCenter

    # Plot the amp positions as a test
    fig, ax = plt.subplots(1)

    # Colors 
    colors = {"RU" : "blue", 
              "RL" : "cyan",
              "LL" : "yellow",
              "LU" : "red"}

    # Plot the polygons
    for amp, rects in six.iteritems(amp_corners):
        for rect in rects:
            corners = []
            for x,y in zip(rect[0], rect[1]):
                corners.append([x, y])
            poly = patches.Polygon(corners, closed=True, alpha=0.25, 
                                   edgecolor="none", facecolor=colors[amp])
            ax.add_patch(poly)

    # Now overplot the fibers
    fibn = []
    x = []
    y = []
    fn = "/afs/ipp-garching.mpg.de/home/d/dfarrow/hetdex/code/virus_config/IFUcen_files/IFUcen_HETDEX.txt"
    with open(fn, 'r') as fp:
        for line in fp:
            if "#" in line:
                continue
            if len(line.split()) < 6:
                continue
            else:
                els = line.strip().split()
                fibn.append(int(els[0]))
                x.append(float(els[1]))
                y.append(float(els[2]))

    fibn = array(fibn)
    x = array(x)
    y = array(y)

    for amp, range_ in six.iteritems(amp_fiber_ranges):
        sel = (fibn >= range_[0]) & (fibn <= range_[1])

        plt.plot(x[sel], y[sel], "o", markerfacecolor=colors[amp],
                 markeredgecolor="k")

        plt.plot(x[sel] - 1.270 , y[sel] + 0.730, "o", markerfacecolor=colors[amp],
                 markeredgecolor="k")

        plt.plot(x[sel] - 1.270 , y[sel] - 0.730, "o", markerfacecolor=colors[amp],
                 markeredgecolor="k")



    plt.ylim(-35, 35)
    plt.xlim(-35, 35)
    plt.xlabel("x IFU (arcseconds)")
    plt.ylabel("y IFU (arcseconds)")
    plt.show()
