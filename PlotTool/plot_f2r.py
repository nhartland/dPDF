# nnpydf - nh 08/14
import sys, os
#nnpydf imports
from thpredictions import ThPredictions
from commondata import CommonData

from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator

from operator import add, sub
import numpy as np
import math

ylims = [0.9,1.1]
xmin = 1E-3
xmax = 0.8
xch  = 0.1

labels = {
  "F2R1":"$c(x)$ (1 GeV$^2$)",
  "F2R10":"$c(x)$ (10 GeV$^2$)",
  "F2R100":"$c(x)$ (100 GeV$^2$)",
  "F2R1000":"$c(x)$ (1000 GeV$^2$)"
}

def compute_MMHT(xvals):
  xp = 0.03
  N  = 0.589
  c1 = -0.116
  c2 = -0.384
  c3 = 0.0489E-8

  yvals = []
  for x in xvals:
    y=0
    if x < xp:
      y = (1+0.01*N)*(1.0 + 0.01*c1*pow(math.log(xp/x),2))
    else:
      y = (1+0.01*N)*(1.0 + 0.01*c2*pow(math.log(x/xp),2) + 0.01*c3*pow(math.log(x/xp),20))
    yvals.append(y)
  return yvals

def cx_plot(prefix, cData, theories):
  gs = gridspec.GridSpec(1, 2)
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
  ax.set_ylim(ylims)
  lax.set_ylim(ylims)
  ax.set_xlim ([xch, xmax])
  lax.set_xlim([xmin, xch])

  plt.setp(ax.get_yticklabels(), visible=False)
  lax.set_xscale('log')
  ax.xaxis.set_major_locator(MaxNLocator(4,prune='lower'))

  ax.set_xlabel("$x$")
  lax.set_xlabel("$x$")
  lax.set_ylabel("$F^2_D / (F^2_n + F^2_p)$")

  for theory in theories:
    CVup = list(map(add, theory.CV, theory.error))
    CVdn = list(map(sub, theory.CV, theory.error))
    ax.plot(cData.kin1, theory.CV, alpha=0.8)
    lax.plot(cData.kin1, theory.CV, alpha=0.8)
    ax.fill_between(cData.kin1, CVdn, CVup, alpha=0.2)
    lax.fill_between(cData.kin1, CVdn, CVup, alpha=0.2, label = labels[cData.setname])

  # Legend
  ax.plot(cData.kin1, compute_MMHT(cData.kin1), alpha=0.8)
  lax.plot(cData.kin1, compute_MMHT(cData.kin1), alpha=0.8, label = "MMHT14 NNLO")
  ax.hlines (1, xch, xmax,  linewidth=1, color = 'k')
  lax.hlines(1, xmin, xch,  linewidth=1, color = 'k')
  legend = lax.legend(fontsize=10, loc='best')
  legend.get_frame().set_alpha(0.7)
  fig.savefig(prefix+cData.setname+'.pdf')


# Fetch datasets
root = sys.argv[1]
datroot = "../data/commondata/"
thrroot  = root + "thr/"

# Make dir for results
plotDir = root + "figures/"
if not os.path.exists(plotDir):
  os.mkdir(plotDir)

for dataset in os.listdir(datroot):
  if dataset[0:6] == "DATA_F" and "_D" not in dataset:
    cDat = CommonData(datroot+dataset)
    prefix = plotDir + cDat.setname + "/"
    c_x  = ThPredictions(thrroot + "TH_" + cDat.setname + ".dat") 
    plot = cx_plot(plotDir, cDat, [c_x])

