"""

Script to plot the sn??.use files 

AUTHOR: Daniel Farrow (MPE)

"""
import matplotlib.pyplot as plt
from matplotlib.colors import TABLEAU_COLORS
from argparse import ArgumentParser
from hetdex_api.flux_limits.flim_models import read_karl_file


parser = ArgumentParser(description="Plot sn.?.?.use files")
parser.add_argument("files", nargs="+", help="Files to plot")
opts = parser.parse_args()


linestyles = ["-", "--", ":", "-."]
first = True
for linestyle, file_ in zip(linestyles, opts.files):

   waves, f50, compl_curves, fluxes = read_karl_file(file_)

   plt.figure(1)
   plt.plot(waves, f50, linestyle=linestyle, label=file_, color="k")
   plt.xlabel("Wavelength (A)")
   plt.ylabel("$10^{-17}$ erg/s/cm$^2$")
   plt.legend()

   plt.figure(2)
   for i, w in enumerate(waves):
       norm = 1.0/max(compl_curves[i, :])
       plt.plot(fluxes, norm*compl_curves[i, :], linestyle=linestyle,
                color=list(TABLEAU_COLORS)[i], label="{:4.0f}".format(w))

   if first:
       plt.legend(frameon = False)
       first = False
   plt.xlabel("Flux [$10^{-17}$ erg/s/cm$^2$]")
   plt.ylabel("Completeness")

plt.show()

