# C Jerome K. Butler July 2001
#
# Make changes that reflect the directory structure below
INCLUDEPATH	= -I$(HOME)/software/bscan/include
LIBPATH		= -L$(HOME)/software/bscan/lib


CC	= gcc
FC	= g77
CFLAGS	= -g -Wall
LFLAGS	= -o Bscan
LIBS	= -lJkb -lcomplx -llinear -lclapack -lblas12f -lF77 -lI77 -lm

HEADERS	= Green.h Dft.h
OBJS	= Bscan.o Green.o Dft.o Znl.o Fields.o
 
Bscan :	$(OBJS) $(HEADERS)
	$(CC) $(LFLAGS) $(OBJS) $(LIBPATH) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDEPATH) -c $*.c

clean :
	rm -f *.o Bscan *.prn out*
