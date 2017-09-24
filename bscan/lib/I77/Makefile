#---------------------------------------------------------------------
#
#  Jerome K. Butler 
#  November 21, 1996
#  Last updated Feb 5, 2001 
#
#######################################################################
CC		= gcc
CFLAGS		= -fPIC -c -O2
#
#---------------------------------------------------------

OBJ =	Version.o backspace.o close.o dfe.o dolio.o due.o endfile.o err.o \
	fmt.o fmtlib.o iio.o ilnw.o inquire.o lread.o lwrite.o open.o \
	rdfmt.o rewind.o rsfe.o rsli.o rsne.o sfe.o sue.o typesize.o uio.o \
	util.o wref.o wrtfmt.o wsfe.o wsle.o wsne.o xwsne.o
#-------------------------------------------------------
library	: $(OBJ)
	gcc -shared -Wl,-soname,libI77.so.1 -o libI77.so.1.0 $(OBJ)

install : 
	mv libI77.so.1.0 ../
	ln -sf libI77.so.1.0 libI77.so.1
	mv libI77.so.1 ../
	ln -sf libI77.so.1 libI77.so
	mv libI77.so ../

clean:
	rm -f *.o

.c.o: 
	$(CC) $(CFLAGS) -c $*.c
