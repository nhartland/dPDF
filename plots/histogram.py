#!/usr/bin/env python3
import sys
from os.path import basename
import os
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import math
import subprocess

colours = ['red', 'blue', 'green']

# Standard Data/Theory comparison
class stdhist:
	def __init__(self, data, labels):
		self.fig = plt.figure()
		self.ax = plt.gca()

		total_data = []
		for set in data:
			for point in set:
				total_data.append(point)
		sorted_data  = sorted(total_data)
		n_datapoints = len(sorted_data)

		# Freedman–Diaconis
		Q1 = sorted_data[int(n_datapoints*0.25)]
		Q3 = sorted_data[int(n_datapoints*0.75)]
		IQR = Q3 - Q1
		h = 2.0*IQR/pow(n_datapoints, 1.0/3.0)
		nbins = math.ceil((max(sorted_data) - min(sorted_data))/h)	

		icol = 0
		for iset in range(0,len(data)):
			self.ax.hist(data[iset], nbins, [min(sorted_data), max(sorted_data)], facecolor=colours[icol], alpha=0.50, label = labels[iset], histtype = 'stepfilled')
			self.ax.legend()
			self.ax.grid(True)
			icol = icol + 1

	def savePlot(self,plotname):
		self.fig.savefig(plotname)


def readfile(filename):
	with open(filename, 'r') as f:
		lines = f.readlines()
	return [float(e.strip()) for e in lines]


# Process set
path = sys.argv[1]
basepath = basename(path)

os.system('rm -rf ./' + basepath)
os.mkdir(basepath) 
os.mkdir(basepath+'/dat/')
os.mkdir(basepath+'/plt/')
basecmd = 'tail -qn1 '+ path + '/erf/replica_*.dat'
os.system(basecmd + '| awk \'{print $1}\' > ./' + basepath + '/dat/TL.dat')
os.system(basecmd + '| awk \'{print $2}\' > ./' + basepath + '/dat/ET.dat')
os.system(basecmd + '| awk \'{print $3}\' > ./' + basepath + '/dat/EV.dat')
os.system(basecmd + '| awk \'{print $4}\' > ./' + basepath + '/dat/C2.dat')

TL = readfile(basepath + '/dat/TL.dat')
ET = readfile(basepath + '/dat/ET.dat')
EV = readfile(basepath + '/dat/EV.dat')
C2 = readfile(basepath + '/dat/C2.dat')

stdhist([TL], ["Training Length"]).savePlot(basepath+"/plt/TL.pdf")
stdhist([C2], ["$\chi^2$"]).savePlot(basepath+"/plt/C2.pdf")
stdhist([ET,EV], ["Training", ["Validation"]]).savePlot(basepath+"/plt/ETEV.pdf")
