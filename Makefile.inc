################################################################################
#
#	dPDF Makefile includes
#
################################################################################

CC= g++
LD= g++

SHELL=/bin/bash

# APPLCOMB paths
RESULTSDIR=   -D  RESULTS_PATH="./results/" -D  DATA_PATH="./data/"
DATABASEDIR=  -D  DB_PATH="./db/"
ALLDIR= $(RESULTSDIR) $(DATABASEDIR)

#LHAPDF
LHAPDFINCS = -I$(shell lhapdf-config --prefix)/include
LHAPDFDIR  = $(shell lhapdf-config --prefix)/lib
LHAPDFLIBS = -L$(LHAPDFDIR) -lLHAPDF

# APFEL
APFELINCS = $(shell apfel-config --cppflags) 
APFELLIBS = $(shell apfel-config --ldflags) 

# libnnpdf
NNPDFINCLUDE= $(shell nnpdf-config --cppflags)
NNPDFLIBS= $(shell nnpdf-config --ldflags)

# libconfig
CONFIGINCLUDE=$(shell pkg-config libconfig --cflags)
CONFIGLIBS=$(shell pkg-config libconfig --libs) -lconfig++ 

# gsl
GSLINCLUDE=$(shell gsl-config --cflags)
GSLLIBS=$(shell gsl-config --libs)

# additional libraries to be included 
PRJLDFLAGS = $(LHAPDFLIBS) $(APPLCLIBS) $(ROOTLIBS) $(APFELLIBS) $(NNPDFLIBS) $(GSLLIBS) $(CONFIGLIBS) $(APFELLIBS) -lsqlite3

# scheduling and optimization options (such as -DSSE -DSSE2 -DP4)
PRJCXXFLAGS = -Wall -O3 -std=c++11 $(ALLDIR) $(LHAPDFINCS) $(NNPDFINCLUDE) $(GSLINCLUDE) $(CONFIGINCLUDE) $(APFELINCS)-I./inc
