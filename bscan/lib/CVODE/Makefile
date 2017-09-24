#---------------------------------------------------------------------
#
#  Jerome K. Butler 
#  November 21, 1996
#  Last updated April 1999 
#
#######################################################################
CC		= gcc
CFLAGS		= -O -fPIC
INCLUDEPATH	= -I/home/jkb/include
#
#---------------------------------------------------------
# The Object Files
#---------------------------------------------------------
OBJS = \
	cvode.o \
	cvdense.o \
	dense.o \
	cvband.o \
	band.o \
	cvdiag.o \
	cvspgmr.o \
	spgmr.o \
	iterativ.o \
	cvbandpre.o \
	nvector.o \
	llnlmath.o
#---------------------------------------------------------
library	: $(OBJS)
	gcc -shared -Wl,-soname,libcvode.so.1 -o libcvode.so.1.0 $(OBJS)

install : 
	ln -sf libcvode.so.1.0 libcvode.so.1
	ln -sf libcvode.so.1 libcvode.so
	mv -f libcvode.so ../
	mv -f libcvode.so.1 ../
	mv -f libcvode.so.1.0 ../

clean:
	rm -f *.o
	rm -f libcvode.so.1.0
	rm -f libcvode.so.1
	rm -f libcvode.so

.c.o: 
	$(CC) $(CFLAGS) $(INCLUDEPATH)  -c $*.c

