include Makefile.inc

CXXFLAGS = 	$(PRJCXXFLAGS) 
LDLIBS   =  $(PRJLDFLAGS)

VPATH=./src
MAIN = dPDF compute_f2 compute_th fkgen
TEST = test1 test2
#ftdy_hcx

.PHONY: all clean
	
all: $(MAIN)
test: $(TEST)

dPDF: dPDF.o filter.o cmaes.o ns_network.o fastaddchi2.o
compute_th: compute_th.o filter.o ns_network.o fastaddchi2.o
compute_f2: compute_f2.o filter.o ns_network.o fastaddchi2.o
fkgen: fkgen.o

test1: test1.o
test2: test2.o ns_network.o

clean:
	-$(RM) -f $(MAIN) $(TEST)
	-$(RM) *.o
