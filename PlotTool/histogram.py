#!/usr/bin/env python3
import sys
from os.path import basename
import os
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import math
import subprocess
from statistics import mean, stdev
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
		nbins = max(10, nbins)
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

def printInfo(label, data):
	data.sort()
	Q14 = data[int(len(data)*0.14)]
	Q86 = data[int(len(data)*0.86)]
	QMd = Q14 + (Q86 - Q14) / 2.0
	print(" ***** " + label + " ***** ")
	print("Arithmetic mean: " + str(mean(data)) + " ± " + str(stdev(data)))
	print("       68% C.I.: " + str(QMd) + " ± " + str(Q86 - QMd))

# Process set
path = sys.argv[1]
if not os.path.exists(path + "/figures"): os.mkdir(path + "/figures") 
basecmd = 'tail -qn1 '+ path + '/erf/replica_*.dat'
os.system(basecmd + '| awk \'{print $1}\' > ' + path + '/erf/TL.dat')
os.system(basecmd + '| awk \'{print $2}\' > ' + path + '/erf/ET.dat')
os.system(basecmd + '| awk \'{print $3}\' > ' + path + '/erf/EV.dat')
os.system(basecmd + '| awk \'{print $4}\' > ' + path + '/erf/C2.dat')

TL = readfile(path + '/erf/TL.dat')
ET = readfile(path + '/erf/ET.dat')
EV = readfile(path + '/erf/EV.dat')
C2 = readfile(path + '/erf/C2.dat')

printInfo("Chi-squared", 		C2)
printInfo("Training Error", 	ET)
printInfo("Validation Error", 	EV)
printInfo("Training Length", 	TL)

stdhist([TL], ["Training Length"]).savePlot(path+"/figures/TL.pdf")
stdhist([C2], ["$\chi^2$"]).savePlot(path+"/figures/C2.pdf")
stdhist([ET,EV], ["Training", ["Validation"]]).savePlot(path+"/figures/ETEV.pdf")
