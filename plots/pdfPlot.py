#! /usr/bin/python
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator

from operator import add, sub
import sys, os, math
import numpy as np
import glob

normalise = True
colours = ['r', 'b', 'g']
pdfnames = ['Singlet', 'Gluon', 'V', 'V8', 'T8']
mxpdf = 3

plotReplicas = True

# Choice of PDF
for pdfidx in xrange(0,mxpdf):

  # Setup gridspec
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

  # Gridlines
  lax.xaxis.grid(True)
  lax.yaxis.grid(True)

  # Colorings
  icol = 0

  # Axis formatting
  plt.setp(ax.get_yticklabels(), visible=False)
  lax.set_xscale('log')
  ax.xaxis.set_major_locator(MaxNLocator(5,prune='lower'))

  lax.set_ylabel("r" + pdfnames[pdfidx])
  ax.set_xlabel("x")
  lax.set_xlabel("x")

  for idat in xrange(1,len(sys.argv)):

    # Fetch path for replica files  
    indir = sys.argv[idat]
    basename = os.path.splitext(sys.argv[idat])[0]

    # Verify paths
    if os.path.exists(indir) == False:
      print "Error: source directory" + indir + " not found!"
      sys.exit()

    xvals = []
    yvals = []

    irep = 0
    nreps = len(os.listdir(indir))
    interval = math.floor((nreps-0.68*nreps)/2)

    print nreps, " Replicas found and processing"
    print interval, " replicas excluded on either side for 68cl"

    for fn in glob.glob(indir+"*.dat"):
      infile = open(fn, 'rb')
      datafile = open(fn, 'rb')

      ix = 0 # x point index
      for line in datafile:
        linesplit = line.split()

        if irep == 0:
          xvals.append(float(linesplit[0]))
          yvals.append([])

        yvals[ix].append(float(linesplit[pdfidx+1]))
        ix = ix + 1

      irep = irep + 1

      # Plot replica
      repyVals = []
      for yListVal in yvals:
        repyVals.append(yListVal[-1])

      if plotReplicas==True:
        ax.plot(xvals,repyVals, color = colours[icol], alpha=0.05)
        lax.plot(xvals,repyVals, color = colours[icol], alpha=0.05)

    # Compute averages
    yCV = []
    yER = []
    y84 = []
    y16 = []
    for replicas in yvals:
      srtreps = np.sort(replicas)
      yCV.append(np.mean(srtreps))
      yER.append(np.std(srtreps))
      y84.append(srtreps[int(nreps-interval-1)])
      y16.append(srtreps[int(interval)])

    CVup = map(add, yCV, yER)
    CVdn = map(sub, yCV, yER)


    # Central values
    ax.plot(xvals,yCV, color = colours[icol], alpha=0.8, label = basename)
    lax.plot(xvals,yCV, color = colours[icol], alpha=0.8, label = basename)

   # error bars
    ax.fill_between(xvals, CVdn, CVup, alpha=0.2,
                          facecolor=colours[icol], linewidth=0, color=colours[icol])
                          
   # error bars
    lax.fill_between(xvals, CVdn, CVup, alpha=0.2,
                          facecolor=colours[icol], linewidth=0, color=colours[icol])

  # 68CL
    ax.plot(xvals,y84, color = colours[icol], alpha=0.8, linestyle='-.')
    lax.plot(xvals,y84, color = colours[icol], alpha=0.8, linestyle='-.')

      # 68CL
    ax.plot(xvals,y16, color = colours[icol], alpha=0.8, linestyle='-.')
    lax.plot(xvals,y16, color = colours[icol], alpha=0.8, label = basename+" 68% C.I", linestyle='-.')

    icol=icol+1

  # set limits
  ax.set_xlim([0.1, 2])
  lax.set_xlim([1E-5,0.1])

  # set limits
  ax.set_ylim([-0.5, 3])
  lax.set_ylim([-0.5,3])

  # Legend
  legend = lax.legend(fontsize=10, loc='best')
  legend.get_frame().set_alpha(0.7)

  fig.savefig('pdfcomp'+pdfnames[pdfidx]+'.pdf')




