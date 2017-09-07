#! /usr/bin/python
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator
import matplotlib.patches as mpatches

from operator import add, sub
import sys, os, math
import numpy as np
import glob

colours = ['b', 'r', 'g']
pdfnames = ['photon', 'Singlet', 'Gluon', 'V', 'V3', 'V8','V15','V24','V35','T3','T8','T15','T24','T35']

# Results directories
pdf_paths  = sys.argv[1:] 
pdf_names = [pdf_path.rsplit('/')[-2] for pdf_path in pdf_paths]

# Make dir for results
plotDir = pdf_paths[0] + "figures/"
if not os.path.exists(plotDir):
  os.mkdir(plotDir)

pdf_figures = []
for ipdf in range(0,len(pdfnames)):
  gs = gridspec.GridSpec(1, 2) #, height_ratios=[1, 1])
  gs.update(wspace=0.00, hspace=0.00)
  w, h = plt.figaspect(0.5)

  # Setup figure
  fig = plt.figure(figsize=(w,h))
  pdf_figures.append(fig)
  lax = fig.add_subplot(gs[0])
  ax  = fig.add_subplot(gs[1])

  # Axis formatting
  ax.set_xlim([0.1, 1])
  lax.set_xlim([1E-5,0.1])
  plt.setp(ax.get_yticklabels(), visible=False)
  lax.set_xscale('log')
  ax.xaxis.set_major_locator(MaxNLocator(5,prune='lower'))
  lax.set_ylabel("r" + pdfnames[ipdf])

  for axis in fig.get_axes():
    axis.xaxis.grid(True)
    axis.yaxis.grid(True)
    axis.set_xlabel("x")
    axis.set_ylim([-1,5])

for iset in range(0, len(pdf_paths)):
  pdfroot = pdf_paths[iset] + 'pdf/'
  print("Reading dataset: " + pdf_names[iset])

  for replica_file in os.listdir(pdfroot):
    replica = np.loadtxt(pdfroot+replica_file)

    for ipdf in range(0,len(pdfnames)):
      for axis in pdf_figures[ipdf].get_axes():
        axis.plot(replica[:,0],replica[:,ipdf+1], color=colours[iset], alpha=0.05)

for ipdf in range(0,len(pdfnames)):
  print("Writing " + pdfnames[ipdf])

  legend_handles = [] # Setup legend
  for iset in range(0, len(pdf_names)):
    legend_handles.append(mpatches.Patch(color=colours[iset], label=pdf_names[iset]))
  pdf_figures[ipdf].get_axes()[0].legend(handles = legend_handles)
  pdf_figures[ipdf].savefig(plotDir + pdfnames[ipdf]+'.pdf')