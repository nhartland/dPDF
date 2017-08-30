# Data plotter in Python
# nh 07/2014

# sys imports
import sys, os

#nnpydf imports
from thpredictions import ThPredictions
from commondata import CommonData
from stdplot import stdplot
from numpy import sqrt

# base class for data plots
class dataplot:
  def __init__(self, data):
    
    if isinstance(data, CommonData) == False:
      print("Error: argument passed to dataplot is not data")
      quit()

    # Copy data
    self.setname = data.setname
    
    self.data = list(data.data)
    self.error = list(data.error)

    self.kin1 = list(data.kin1)
    self.kin2 = list(data.kin2)
    self.kin3 = list(data.kin3)
    
    self.pdf = []
    self.theory = []
    self.theoryerror = []

  def addTheory(self, theory):
      if isinstance(theory, ThPredictions) == False:
        print("Error: argument passed to dataplot.addTheory is not theory")
        quit()
      
      if theory.setname != self.setname:
        print("Error: theory and data setnames do not match!",theory.setname,self.setname)
        quit()

      # Copy data
      self.pdf.append(theory.pdfset)
      
      self.theory.append(list(theory.CV))
      self.theoryerror.append(list(theory.error))

# *********************************************************************************************************************


# Data plotting with no kinematics info
class dataplot_nokin(dataplot):
  def __init__(self, data):
    dataplot.__init__(self, data)
    
    # Initialise plot
    xaxis = range(1,len(self.data)+1)
    self.stdplot = stdplot(xaxis, self.data, self.error, data.setname, xlab = "datapoint")

  def addTheory(self,theory):
    dataplot.addTheory(self,theory)
    
    # theory values
    CV = self.theory[len(self.theory)-1]
    ER = self.theoryerror[len(self.theory)-1]
    
    self.stdplot.addTheory(CV,ER, theory.pdfset)
  
  def export(self, prefix = "./"):
    self.stdplot.savePlot(prefix+ self.setname)


# *********************************************************************************************************************


# Data plotting with kinematics dependance - base class
class dataplot_kin_base(dataplot):
  def __init__(self, data, kinx, kiny = 4, kinz = 4, sqrtx = False):
    dataplot.__init__(self, data)
    
    # Set chosen binning
    self.kinx = kinx
    self.kiny = kiny
    self.kinz = kinz
    
    # Square root the x-axis
    self.sqrtx = sqrtx
    
    # Kinematics dictionary
    self.kinDict = [None, self.kin1, self.kin2, self.kin3, [None for x in range(len(self.data))]]
    
    # Get y kinematics
    self.kinyVals = list(set(self.kinDict[self.kiny]))
    self.kinyVals.sort(reverse=True) # Reverse kin for DIS
    
    # Get z kinematics
    self.kinzVals = list(set(self.kinDict[self.kinz]))
    self.kinzVals.sort()
    
    # Number of bins
    self.ny = len(self.kinyVals)
    self.nz = len(self.kinzVals)
    
    # Plan and index lists - default values
    self.pltindx = [[[None,None] for x in range(self.ny)] for x in range(self.nz)]
    self.dataidx = [[None,None] for x in range(len(self.data))]
    
    # Split data arrays
    self.splitData = [[ None for x in range(self.ny)] for x in range(self.nz)]
    self.splitErr = [[ None for x in range(self.ny)] for x in range(self.nz)]
    self.splitKin = [[ None for x in range(self.ny)] for x in range(self.nz)]
    
    # Split theory arrays
    self.splitTheory = [[ [] for x in range(self.ny)] for x in range(self.nz)]
    self.splitTheoryErr = [[ [] for x in range(self.ny)] for x in range(self.nz)]
    self.splitTheoryKin = [[ [] for x in range(self.ny)] for x in range(self.nz)]
    
    # Loop over indices
    for zIdx in range(self.nz):
      for yIdx in range(self.ny):
        plotkiny = self.kinyVals[yIdx]
        plotkinz = self.kinzVals[zIdx]
        
        # set index
        self.pltindx[zIdx][yIdx] = [plotkinz,plotkiny]
        
        # Subplot data
        subDat = []
        subKin = []
        subErr = []
        
        #Loop through datapoints
        for i in range(0,len(self.data)):
          datkiny = self.kinDict[self.kiny][i]
          datkinz = self.kinDict[self.kinz][i]
          
          #If data kinematics match
          if datkiny == plotkiny and datkinz == plotkinz:
            self.dataidx[i] = [zIdx,yIdx]
            subDat.append(self.data[i])
            subErr.append(self.error[i])
            
            if sqrtx == True:
              subKin.append(sqrt(self.kinDict[self.kinx][i]))
            else:
              subKin.append(self.kinDict[self.kinx][i])
      
        # Push into data array
        if len(subDat) != 0:
          self.splitData[zIdx][yIdx] = subDat
          self.splitErr[zIdx][yIdx] = subErr
          self.splitKin[zIdx][yIdx] = subKin
  
  
  def addTheory(self,theory):
    dataplot.addTheory(self,theory)
    
    # theory values
    CV = self.theory[len(self.theory)-1]
    ER = self.theoryerror[len(self.theory)-1]
    
    # Loop over indices
    for zIdx in range(self.nz):
      for yIdx in range(self.ny):
        
        # Subplot data
        subDat = []
        subKin = []
        subErr = []
        
        #Loop through datapoints
        for i in range(0,len(CV)):
          if self.dataidx[i] == [zIdx, yIdx]:
            subDat.append(CV[i])
            subErr.append(ER[i])
            
            if self.sqrtx == True:
              subKin.append(sqrt(self.kinDict[self.kinx][i]))
            else:
              subKin.append(self.kinDict[self.kinx][i])
        
        # Push into theory array
        if len(subDat) != 0:
          self.splitTheory[zIdx][yIdx].append(subDat)
          self.splitTheoryErr[zIdx][yIdx].append(subErr)
          self.splitTheoryKin[zIdx][yIdx].append(subKin)

# *********************************************************************************************************************

# Data plotting with kinematics dependance
class dataplot_kin(dataplot_kin_base):
  def __init__(self, data, kinx, kiny = 4, kinz = 4, logx = False, logy = False, sqrtx = False,
               xlab = "x", ylab="y", zlab = "z"):
    dataplot_kin_base.__init__(self, data, kinx, kiny, kinz, sqrtx)
    
    # Plan and index lists - default values
    self.stdplots = [[0 for x in range(self.ny)] for x in range(self.nz)]
    
    # Loop over y indices
    for zIdx in range(self.nz):
      for yIdx in range(self.ny):
        if self.splitData[zIdx][yIdx] != None:
          plotkiny = self.kinyVals[yIdx]
          plotkinz = self.kinzVals[zIdx]
          
          # Initialise plot
          self.stdplots[zIdx][yIdx] = stdplot(self.splitKin[zIdx][yIdx], self.splitData[zIdx][yIdx], self.splitErr[zIdx][yIdx], data.setname, logx, logy, xlab )
          # Plot Text
          kintitle = None
          
          if plotkiny != None:
            kintitle = ylab+' = ' + str(round(plotkiny,4))
          if plotkinz != None:
            kintitle = kintitle +', '+zlab+' = '+str(round(plotkinz,2))
          
          if kintitle != None:
            self.stdplots[zIdx][yIdx].fig.text(0.90,0.92,kintitle, fontsize=12, horizontalalignment='right')
        else:
          self.stdplots[zIdx][yIdx] = None
  
  
  def addTheory(self,theory):
    dataplot_kin_base.addTheory(self,theory)

    for yIdx in range(0, self.ny):
      for zIdx in range(0, self.nz):
        if self.stdplots[zIdx][yIdx] != None:
          self.stdplots[zIdx][yIdx].addTheory(self.splitTheory[zIdx][yIdx][-1],self.splitTheoryErr[zIdx][yIdx][-1],theory.pdfset )
  
  def export(self, prefix = "./"):
    for yIdx in range(0, self.ny):
      for zIdx in range(0, self.nz):
        if self.stdplots[zIdx][yIdx] != None:
          self.stdplots[zIdx][yIdx].savePlot(prefix+self.setname+"_"+str(yIdx)+"_"+str(zIdx))