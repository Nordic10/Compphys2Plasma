FC	=	gfortran
#LFLAGS = 	-O -L$(COMPHY)/lib -lnumer
LFLAGS 	= 	-g -Llib -lnumer
FFLAGS 	= 	-c -g

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs) 
CXXFLAGS  += $(ROOTCFLAGS)
GLIBS      = $(ROOTGLIBS)

ROOTFLAGS   = $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS) 
#P5640FLAGS  = -L${P5640LIB}/lib -lP5640  -I${P5640LIB}
GSLFLAGS     = -lgsl -lgslcblas




all: Vgen VgenVec


Vgen: Vgen.C
	g++ -g -Wall -oVgen Vgen.C $(ROOTFLAGS) $(GSLFLAGS)

VgenVec: VgenVec.C
	g++ -g -Wall -oVgenVec VgenVec.C $(ROOTFLAGS) $(GSLFLAGS)

clean:
	rm -f VgenVec *.o *.so *.pcm *.d *~
