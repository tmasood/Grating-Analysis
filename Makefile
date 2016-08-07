# Copyright 2015 Quantum Designs LLC, Taha Masood, Johannes Tausch
# and Jerome Butler
#
# Make changes that reflect the directory structure below
INCLUDEPATH	= -I/usr/include
LIBPATH		= -L/usr/lib

CC	= g++
FC	= g++
CFLAGS	= -g -Wall -fno-common
LFLAGS	= -o gratanal
LIBS	= -llapack -lm -lblas -lQP -lbessel -lgfortran

HEADERS	= layer.h structure.h grating.h ssystem.h newton.h gtoothpnl.h gtoothdmn.h
OBJS	= gratanal.o structure.o layer.o readfile.o grating.o ssystem.o readlayer.o initialize.o newtoncharmatrix.o znewton.o gtoothpnl.o defpnl.o gtoothdmn.o defdmn.o getcollocpts.o indxpnl.o refinepnls.o gettannrm.o initcalcp.o calcp.o calcrhs.o calcrhs1.o getlogweights.o calcfc.o calcfc1.o calcmoments.o initsheets.o deconvolve.o convolve.o setuptransops.o translated2n.o translated2nr.o subext.o infnorm.o nextgam.o dumpnls.o getnullspace.o translatedvw.o calcsolnpt.o printsolution.o calcrhsint.o domainnr.o calcsolint.o printsolint.o testsol.o

gratanal:	$(OBJS) $(HEADERS)
	$(FC) $(LFLAGS) $(OBJS) $(LIBPATH) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDEPATH) -c $<

clean :
	rm -f core *.o gratanal *.prn out*

install:
	mv gratanal /home/tmasood/Desktop/modeling/tools/grating/bin
	chmod 755 /home/tmasood/Desktop/modeling/tools/grating/Bscan
