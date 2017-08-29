#!/usr/bin/env python3
import sys
from os.path import basename
import os
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import math
import subprocess
import numpy as np

colours = ['red', 'blue', 'green']
sets = ['BCDMSD', 'SLACD', 'DYE886R', 'NMCPD']
# sets = ['DYE886R']

# Process set
path = sys.argv[1]
basepath = basename(path)

for set in sets:
	# Deuteron predictions
	tpath = './' + basepath + '/dat/TH_'+set+'.dat'
	os.system('paste '+ path + '/thd/'+set+'_replica_*.dat > ' + tpath)
	deuteron = np.genfromtxt(tpath)

	tpath = './' + basepath + '/dat/PR_'+set+'.dat'
	os.system('paste '+ path + '/thp/'+set+'_replica_*.dat > ' + tpath)
	proton = np.genfromtxt(tpath)

	data   = np.genfromtxt(path + '/dat/' + set +'.dat')

	fig, ax = plt.subplots(1)
	ax.errorbar(range(0, len(data[:,0])), data[:,0], yerr = data[:,1], fmt='.', color='black')
	for i in range(deuteron.shape[1]):
	    ax.plot(deuteron[:,i], '-b', alpha=0.1)
	for i in range(proton.shape[1]):
	    ax.plot(proton[:,i], '-r', alpha=0.1)
	plt.show()