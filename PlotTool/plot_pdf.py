#! /usr/bin/python
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator

from operator import add, sub
import sys, os, math
import numpy as np
import glob

colours = ['r', 'b', 'g']
pdfnames = ['photon', 'Singlet', 'Gluon', 'V', 'V3', 'V8','V15','V24','V35','T3','T8','T15','T24','T35']

# Fetch datasets
root = sys.argv[1]
pdfroot = root + "pdf/"

# Make dir for results
plotDir = root + "figures/"
if not os.path.exists(plotDir):
  os.mkdir(plotDir)

replicas = []
for replica_file in os.listdir(pdfroot):
  replica = np.loadtxt(pdfroot+replica_file)
  replicas.append(replica)

pdf_figures = []
for ipdf in range(0,len(pdfnames)):
  gs = gridspec.GridSpec(1, 2) #, height_ratios=[1, 1])
  gs.update(wspace=0.00, hspace=0.00)
  w, h = plt.figaspect(0.5)

  # Setup figure
  fig = plt.figure(figsize=(w,h))
  pdf_figures.append(fig)

  # Axis formatting
  plt.setp(ax.get_yticklabels(), visible=False)
  lax.set_xscale('log')
  ax.xaxis.set_major_locator(MaxNLocator(5,prune='lower'))
  lax.set_ylabel("r" + pdfnames[ipdf])

  lax = fig.add_subplot(gs[0])
  ax  = fig.add_subplot(gs[1])
  ax.set_xlim([0.1, 1])
  lax.set_xlim([1E-5,0.1])

  for axis in fig.get_axes():
    axis.xaxis.grid(True)
    axis.yaxis.grid(True)
    axis.set_xlabel("x")
    axis.set_ylim([-1,5])

for ipdf in range(0,len(pdfnames)):
  for ipdf in range(0,len(pdfnames)):
    for axis in pdf_figures[ipdf].get_axes():
      axis.plot(replica[:,0],replica[:,ipdf+1], color='blue', alpha=1/len(replicas))
  fig.savefig(plotDir + pdfnames[ipdf]+'.pdf')


   # Legend
#   legend = lax.legend(fontsize=10, loc='best')
#   legend.get_frame().set_alpha(0.7)


