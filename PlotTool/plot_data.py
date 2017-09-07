#!/usr/bin/python
# nnpydf - nh 08/14
import os, shutil
#nnpydf imports
from thpredictions import ThPredictions
from commondata import CommonData
from datplot import *


#############################################################################################################################

def genPlotPage(prefix, dataset, cData, theories):

  print("Generating plots for: " + dataset)

  # Make folder
  if not os.path.exists(prefix):
    os.mkdir(prefix)  
  dataPlots = []
  
  # No kinematics, and DIS plots
  dataPlots.append(dataplot_nokin(cData))
  
  # Individual kinematics bins
  if cData.setname.startswith("F2R"):
    dataPlots.append(dataplot_kin(cData,1,2, xlab="x", ylab="$Q^2$", logx=True))
  elif cData.proc.startswith("DIS"):
    dataPlots.append(dataplot_kin(cData,2,1,3, xlab="$Q^2$", ylab="x", zlab="$E_{CM}$"))
  elif cData.setname.startswith("DYE886R"):
    dataPlots.append(dataplot_kin(cData,1, xlab="rapidity"))
  
  # # Add theory predictions
  for theory in theories:
    for plot in dataPlots:
      plot.addTheory(theory)
  
  for plot in dataPlots:
    plot.export(prefix)

#################################################################################################################################

# Results directories
pdf_paths  = sys.argv[1:] 
pdf_names = [pdf_path.rsplit('/')[-2] for pdf_path in pdf_paths]

# Fetch datasets
root = sys.argv[1]
datroot = root + "dat/"

# Make dir for results
plotDir = root + "figures/"
if not os.path.exists(plotDir):
  os.mkdir(plotDir)

cData = []
for dataset in os.listdir(datroot):
  if dataset[0:4] == "DATA":
    cData.append(CommonData(datroot+dataset))

for cDat in cData:
  prefix = plotDir + cDat.setname + "/"
  theory = []
  for ipdf in range(0,len(pdf_paths)):
    source = pdf_paths[ipdf] + "thd/TH_" + cDat.setname + ".dat"
    if (os.path.exists(source)):
      theory.append(ThPredictions(source))
      theory[-1].pdfset = pdf_names[ipdf]
  genPlotPage(prefix, cDat.setname, cDat, theory)
