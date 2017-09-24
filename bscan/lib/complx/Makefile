#---------------------------------------------------------------------
#
#  Jerome K. Butler 
#  November 21, 1996
#  Last updated August 4, 2000 
#
#######################################################################
# Change the necessary directory structure.  $HOME is the user
# home directory.

INCLUDEPATH	= $(HOME)/software/bscan/include
LIBDIR		= /usr/lib

CC		= gcc
CFLAGS		= -fPIC -I$(INCLUDEPATH) -O2
#
#---------------------------------------------------------
OBJS =	btan.o \
	cabs1.o \
	cadd.o \
	casin.o \
	ccos.o \
	cdd.o \
	cdiv.o \
	ceq.o \
	cfm.o \
	cfadd.o \
	cfsub.o \
	cguesses.o \
	croot.o \
	csin.o \
	csq.o \
	cmul.o \
	cneg.o \
	csub.o \
	csqrt.o \
	powi.o \
	czeta.o \
	dcd.o \
	fmul.o \
	fcm.o \
	fcadd.o \
	fcsub.o \
	tdiv.o \
	fdiv.o \
	zinv.o \
	dconjg.o \
	ftoc.o \
	simp.o \
	tmul.o \
	zcheb.o \
	zconj.o \
	zsqrt.o \
	zlog.o \
	zexp.o \
	ztan.o \
	zabs.o \
	zabs2.o \
	zdiv.o \
	zromb.o \
	zch.o \
	zsh.o \
	zunounp1.o \
	unit.o \
	zvm.o
#-------------------------------------------------------
library	: $(OBJS)
	$(CC) -shared -Wl,-soname,libcomplx.so.1 -o libcomplx.so.1.0 $(OBJS)

install : 
	mv libcomplx.so.1.0 ${LIBDIR}
	ln -sf libcomplx.so.1.0 libcomplx.so.1
	mv libcomplx.so.1 ${LIBDIR}
	ln -sf libcomplx.so.1 libcomplx.so
	mv libcomplx.so ${LIBDIR}

clean:
	rm -f *.o lib*

.c.o: 
	$(CC) $(CFLAGS) -c $*.c

