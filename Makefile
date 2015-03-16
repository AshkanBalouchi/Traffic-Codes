# special variables
CC = g++ 
##CC     = g++ -Wall -ansi -pedantic -pg 
CFLAGS  = -pg -O3
LFLAGS  = -pg
##CFLAGS  = -O3
##LFLAGS  = 
LIBRARIES = -lm -lgsl -lgslcblas 
#LIBRARIES = -lm -lgsl -lgslcblas -framework -lgfortran -lfftw3
#LIBRARIES = -framework veclib -lm
#cd LIBRARIES = -lfftw3 -lm
#LIBRARIES = -lmpi -lscs_mp

# other variables
##MAINFILE = analyse2 
##MAINFILE = sim2 
##MAINFILE = bitFourier2
##MAINFILE = FourierAnalyse 
##MAINFILE = pattern
##MAINFILE = CheckDDC
##MAINFILE = DDC
##MAINFILE = FindGap1
##MAINFILE = Analyitical2
##MAINFILE = Analytical2nd
##MAINFILE = Analytical3rd
##MAINFILE = AnalyticalFull2
##MAINFILE = FindgapFourrier
##MAINFILE = FindGapEvo
##MAINFILE = bashtest
##MAINFILE = FullAnalyse
##MAINFILE = OpenBoundary
##MAINFILE = FullAnalyseBreak
##MAINFILE = FindGapDis
##MAINFILE = Susceptibility
##MAINFILE = NNBreakCorrelation
##MAINFILE = Chi4
##MAINFILE = GapDis
##MAINFILE = Ising4
##MAINFILE = IsingEnergy
##MAINFILE = IsingEnergySpin
##MAINFILE = IsingEnergyMagnetization
##MAINFILE = VelGap
##MAINFILE = JamTimeEvo
##MAINFILE = JamHist
MAINFILE = JamHistFreeCarsDensity





# OTHERS   = traffic_obj
 OBJS = $(MAINFILE:%=%.o) $(OTHERS:%=%.o) 

# pattern rules to define implicit rule
%.o: %.cc
	$(CC) -c $(CFLAGS) $< 
       
%.o: %.cpp
	g++ -c $(CFLAGS) $<
        
%.o: %.f
	f77 -c $<


# default rule(s)
all: $(MAINFILE)
$(MAINFILE): $(OBJS)
	$(LINK.c) $(LFLAGS) -o exec $^ $(LIBRARIES)

# overridden implict rules (i.e. extra dependencies)
main.o: main.cc  traffic_obj.h Makefile
#qobjects.o: qobjects.cc qobjects.h Makefile
#calc_objective.o: main.cc main.h qobjects.h Makefile
#lapack_interface.o: lapack_interface.h Makefile

        
#depend:
#	mkdep -f depend $(OFLAGS) $(MAINFILE:%=%.cc) $(OTHERS:%=%.cc)
#.PHONY: depend
#include depend



#!/bin/bash
