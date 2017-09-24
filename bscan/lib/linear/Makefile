#---------------------------------------------------------------------
#
#  Jerome K. Butler 
#  November 21, 1996
#  Last updated April 30, 1999 
#
#######################################################################
CC		= gcc
INCLUDEPATH	= /home/taha/software/bscan/include
LIBPATH         = ${HOME}/lib
CFLAGS		= -fPIC -c -I$(INCLUDEPATH) -O
#---------------------------------------------------------
#  The Object Files
#---------------------------------------------------------
OBJS =  cout_.o \
	indexx.o \
	iout_.o \
	lud_.o \
	rout_.o \
	matrixmodify.o \
	vectormodify.o
#----------------------------------------------------------
library	: $(OBJS)
	gcc -shared -Wl,-soname,liblinear.so.1 -o liblinear.so.1.0 $(OBJS)

objects : $(SRC)
	$(CC) $(CFLAGS) $(SRC)

install :
	ln -sf liblinear.so.1.0 liblinear.so.1
	ln -sf liblinear.so.1 liblinear.so
	mv -f liblinear.so ${LIBPATH}
	mv -f liblinear.so.1 ${LIBPATH}
	mv -f liblinear.so.1.0 ${LIBPATH}

clean:
	rm -f *.o
	rm -f liblinear.so.1.0
	rm -f liblinear.so.1
	rm -f liblinear.so

.c.o: 
	$(CC) $(CFLAGS) -I$(INCLUDEPATH)  -c $*.c


