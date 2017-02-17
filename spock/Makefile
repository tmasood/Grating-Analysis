# Copyright 2015 Quantum Designs LLC, Taha Masood, Johannes Tausch
# and Jerome Butler
#
# Make changes that reflect the directory structure below
INCLUDEPATH	= -I/usr/include
LIBPATH		= -L/usr/lib

CC	= gcc
FC	= gcc
CFLAGS	= -g -Wall -fno-common
LFLAGS	= -o spock
LIBS	= -lm

HEADERS	= case.h gamout.h local_complex.h output.h util.h charmatrix.h layer.h modcon.h struct.h
OBJS	= readinputfile.o readcase.o readoutput.o readgamout.o readlayers.o readgradedlayer.o commandoptions.o guideparams.o calcindex.o charmatrix.o search.o czerom.o initialize.o local_complex.o eigeqf.o laymat.o powers.o int4pt.o pfail.o nrmlzf.o width.o ftfeld.o field.o readloopz.o numzlp.o zincr.o commands/libcommands.o matsys/libmatsys.o spock.o

spock:	$(OBJS) $(HEADERS)
	$(FC) $(LFLAGS) $(OBJS) $(LIBPATH) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDEPATH) -c $<

clean :
	rm -f core *.o gratanal *.prn out*

install:
	mv gratanal /home/tmasood/Desktop/modeling/tools/grating/bin
	chmod 755 /home/tmasood/Desktop/modeling/tools/grating/Bscan
