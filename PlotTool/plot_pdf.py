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

for ipdf in range(0,len(pdfnames)):
  gs = gridspec.GridSpec(1, 2) #, height_ratios=[1, 1])
  gs.update(wspace=0.00, hspace=0.00)
  w, h = plt.figaspect(0.5)

  # Setup figure
  fig = plt.figure(figsize=(w,h))
  lax = fig.add_subplot(gs[0])
  ax  = fig.add_subplot(gs[1])

  ax.xaxis.grid(True)
  ax.yaxis.grid(True)
  lax.xaxis.grid(True)
  lax.yaxis.grid(True)

  # Axis formatting
  plt.setp(ax.get_yticklabels(), visible=False)
  lax.set_xscale('log')
  ax.xaxis.set_major_locator(MaxNLocator(5,prune='lower'))
  lax.set_ylabel("r" + pdfnames[ipdf])
  ax.set_xlabel("x")
  lax.set_xlabel("x")

  for replica in replicas:
    ax.plot(replica[:,0],replica[:,ipdf+1], color='blue', alpha=1/len(replicas))
    lax.plot(replica[:,0],replica[:,ipdf+1], color='blue', alpha=1/len(replicas))

  # set limits
  ax.set_xlim([0.1, 1])
  lax.set_xlim([1E-5,0.1])

  # # set limits
  # ax.set_ylim([-0.5, 3])
  # lax.set_ylim([-0.5,3])

#   # Legend
#   legend = lax.legend(fontsize=10, loc='best')
#   legend.get_frame().set_alpha(0.7)

  fig.savefig(plotDir + pdfnames[ipdf]+'.pdf')




