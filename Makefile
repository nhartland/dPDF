include Makefile.inc

CXXFLAGS = 	$(PRJCXXFLAGS) 
LDLIBS   =  $(PRJLDFLAGS)

VPATH=./src
MAIN = dPDF test plotter
#ftdy_hcx

.PHONY: all clean
	
all: $(MAIN)

dPDF: dPDF.o filter.o cmaes.o ns_network.o fastaddchi2.o
test: test.o
plotter: plotter.o filter.o ns_network.o fastaddchi2.o

clean:
	-$(RM) -f $(MAIN)
	-$(RM) *.o
