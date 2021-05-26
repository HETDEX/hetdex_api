"""

Fit a polynomial of S/N cut to
noise multiplier using a file of the
form

#name  s/n ccut0 ccut1  f50x
sn4.8.use 4.8 0.000 0.0000 4.49 4.607
sn5.0.use 5.0 0.000 0.0000 4.60 4.783
sn5.5.use 5.5 0.000 0.0000 4.82 5.282
sn6.0.use 6.0 0.000 0.0000 5.18 5.805
sn6.5.use 6.5 0.000 0.0000 5.61 6.362
sn7.0.use 7.0 0.000 0.0000 6.04 6.963

"""
import matplotlib.pyplot as plt
from numpy import polyfit, polyval, array


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

