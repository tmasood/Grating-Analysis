#---------------------------------------------------------------------
#
#  Jerome K. Butler 
#  November 21, 1996
#  Last updated April 1999 
#
#######################################################################
CC		= gcc
INCLUDEPATH     = /home/taha/software/bscan/include
LIBPATH         = ${HOME}/lib
CFLAGS          = -fPIC -c -I$(INCLUDEPATH) -O
#
#---------------------------------------------------------
# The Object Files
#---------------------------------------------------------
OBJS = \
	Dcheb.o \
	Floc.o \
	Powi.o \
	Simp.o \
	dunounp1.o
#---------------------------------------------------------
library	: $(OBJS)
	gcc -shared -Wl,-soname,libJkb.so.1 -o libJkb.so.1.0 $(OBJS)

install : 
	ln -sf libJkb.so.1.0 libJkb.so.1
	ln -sf libJkb.so.1 libJkb.so
	mv -f libJkb.so ../
	mv -f libJkb.so.1 ../
	mv -f libJkb.so.1.0 ../

clean:
	rm -f *.o
	rm -f libJkb.so.1.0
	rm -f libJkb.so.1
	rm -f libJkb.so

.C.o: 
	$(CC) $(CFLAGS) $(INCLUDEPATH)  -c $<

