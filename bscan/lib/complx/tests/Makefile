CC		= gcc 
INCLUDEPATH	= -I/home/jkb/include
CFLAGS		= -c -g -W -Wall $(INCLUDEPATH)
LFLAGS		= -o Testcasin
LIBPATH		= -L/home/jkb/lib
LIBS		= -lcomplx -lm

Testcasin :	Testcasin.o 
	${CC} ${LFLAGS} Testcasin.o \
	$(LIBPATH) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) -c $<

