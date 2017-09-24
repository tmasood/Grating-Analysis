CC		= g++
INCLUDEPATH	= /home/taha/software/bscan/include
CFLAGS		= -c -g -W -Wall -I$(INCLUDEPATH)
#CFLAGS		= -c -O2 -I$(INCLUDEPATH)
LFLAGS		= -o floctest
LIBPATH		= /home/taha/software/bscan/lib
LIBS		=  -lJkb 

floctest : floctest.o 
	${CC} ${LFLAGS} floctest.o \
	-L$(LIBPATH) $(LIBS)

.C.o:
	$(CC) $(CFLAGS) -c $<

clean :
	rm -f *.o floctest

