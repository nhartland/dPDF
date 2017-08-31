#!/usr/bin/python
# nnpydf - nh 08/14
import sys, os, shutil, glob
#nnpydf imports
from thpredictions import ThPredictions
from commondata import CommonData

from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator

from operator import add, sub
import numpy as np


def cx_plot(prefix, cData, theories):
  gs = gridspec.GridSpec(1, 2) #, height_ratios=[1, 1])
  gs.update(wspace=0.00, hspace=0.00)

  # Aspect ratio
  w, h = plt.figaspect(0.5)

  # Setup figure
  fig = plt.figure(figsize=(w,h))
  lax = fig.add_subplot(gs[0])
  ax  = fig.add_subplot(gs[1])

  # Gridlines
  ax.xaxis.grid(True)
  ax.yaxis.grid(True)
  lax.xaxis.grid(True)
  lax.yaxis.grid(True)

  # Set limits
  ax.set_ylim([0.5,1.25])
  lax.set_ylim([0.5,1.25])
  ax.set_xlim([0.1, 0.8])
  lax.set_xlim([1E-3,0.1])

  plt.setp(ax.get_yticklabels(), visible=False)
  lax.set_xscale('log')
  ax.xaxis.set_major_locator(MaxNLocator(5,prune='lower'))

  ax.set_xlabel("x")
  lax.set_xlabel("x")

  for theory in theories:
    CVup = list(map(add, theory.CV, theory.error))
    CVdn = list(map(sub, theory.CV, theory.error))
    ax.plot(cData.kin1, theory.CV, alpha=0.8)
    lax.plot(cData.kin1, theory.CV, alpha=0.8)
    ax.fill_between(cData.kin1, CVdn, CVup, alpha=0.2)#,facecolor=colours[icol], linewidth=0, color=colours[icol])
    lax.fill_between(cData.kin1, CVdn, CVup, alpha=0.2)#,facecolor=colours[icol], linewidth=0, color=colours[icol])


  fig.savefig(prefix+cData.setname+'.pdf')


# Fetch datasets
root = sys.argv[1]
datroot = "../data/commondata/"
thrroot  = root + "thr/"

# Make dir for results
plotDir = root + "figures/"
if not os.path.exists(plotDir):
  os.mkdir(plotDir)

cData = []
for dataset in os.listdir(datroot):
  if dataset[0:6] == "DATA_F":
    cData.append(CommonData(datroot+dataset))

for cDat in cData:
  prefix = plotDir + cDat.setname + "/"
  c_x  = ThPredictions(thrroot + "TH_" + cDat.setname + ".dat") 
  plot = cx_plot(plotDir, cDat, [c_x])

#     for fn in glob.glob(indir+"*.dat"):
#       infile = open(fn, 'rb')
#       datafile = open(fn, 'rb')

#       ix = 0 # x point index
#       for line in datafile:
#         linesplit = line.split()

#         if irep == 0:
#           xvals.append(float(linesplit[0]))
#           yvals.append([])

#         yvals[ix].append(float(linesplit[pdfidx+1]))
#         ix = ix + 1

#       irep = irep + 1

#       # Plot replica
#       repyVals = []
#       for yListVal in yvals:
#         repyVals.append(yListVal[-1])

#       if plotReplicas==True:
#         ax.plot(xvals,repyVals, color = colours[icol], alpha=0.05)
#         lax.plot(xvals,repyVals, color = colours[icol], alpha=0.05)

#     # Compute averages
#     yCV = []
#     yER = []
#     y84 = []
#     y16 = []
#     for replicas in yvals:
#       srtreps = np.sort(replicas)
#       yCV.append(np.mean(srtreps))
#       yER.append(np.std(srtreps))
#       y84.append(srtreps[int(nreps-interval-1)])
#       y16.append(srtreps[int(interval)])

#     CVup = map(add, yCV, yER)
#     CVdn = map(sub, yCV, yER)


#     # Central values
#     ax.plot(xvals,yCV, color = colours[icol], alpha=0.8, label = basename)
#     lax.plot(xvals,yCV, color = colours[icol], alpha=0.8, label = basename)

#   # 68CL
#     ax.plot(xvals,y84, color = colours[icol], alpha=0.8, linestyle='-.')
#     lax.plot(xvals,y84, color = colours[icol], alpha=0.8, linestyle='-.')

#       # 68CL
#     ax.plot(xvals,y16, color = colours[icol], alpha=0.8, linestyle='-.')
#     lax.plot(xvals,y16, color = colours[icol], alpha=0.8, label = basename+" 68% C.I", linestyle='-.')

#     icol=icol+1

#   # set limits
#   ax.set_xlim([0.1, 2])
#   lax.set_xlim([1E-5,0.1])

#   # set limits
#   ax.set_ylim([-0.5, 3])
#   lax.set_ylim([-0.5,3])

#   # Legend
#   legend = lax.legend(fontsize=10, loc='best')
#   legend.get_frame().set_alpha(0.7)

#   fig.savefig('pdfcomp'+pdfnames[pdfidx]+'.pdf')




