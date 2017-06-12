include Makefile.inc

CXXFLAGS = 	$(PRJCXXFLAGS) 
LDLIBS   =  $(PRJLDFLAGS)

VPATH=./src
MAIN = dPDF test
DEV = appl_optgrid
#ftdy_hcx

.PHONY: all clean
	
all: $(MAIN)

dPDF: dPDF.o filter.o cmaes.o ns_network.o fastaddchi2.o proton.o
test: test.o filter.o cmaes.o ns_network.o fastaddchi2.o proton.o
clean:
	-$(RM) -f $(MAIN)
	-$(RM) *.o
