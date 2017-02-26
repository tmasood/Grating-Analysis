CC=gcc
CFLAGS=-Wall -g -pipe -pedantic
INC =   -I /usr/include
LIBS =  -L /home/tmasood/
LDFLAGS = -lg2c -lm 

OBJS= tst.o gamma.o r9lgmc.o gamlim.o csevl.o inits.o \
alngam.o asyjy.o jairy.o cshch.o ckscl.o cs1s2.o cacai.o \
cairy.o cuni1.o cuni2.o cbuni.o cunhj.o cuchk.o crati.o \
cbknu.o cunik.o casyi.o cmlri.o cseri.o cuoik.o cwrsk.o \
gamln.o xgetua.o xersve.o j4save.o xerhlt.o xercnt.o \
fdump.o xerprn.o xermsg.o i1mach.o r1mach.o besy0.o \
besy1.o besynu.o yairy.o besj0.o besj1.o \
cbesj.o besy.o besyp.o besi0e.o besi1e.o besi0.o besi1.o \
besi.o besip.o besknu.o besk0e.o besk1e.o asyik.o besk0.o \
besk1.o besk.o beskp.o matmul.o matmul2.o \
caldet.o cbinu.o matsolv.o sgetrf.o \
sgetri.o xerbla.o ilaenv.o sgemv.o sgemm.o strsm.o sswap.o \
strtri.o strmm.o lsame.o strti2.o sgetf2.o slaswp.o isamax.o \
sger.o sscal.o strmv.o ctransfer.o assign_region.o \
c_compareg.o c_comparege.o c_comparel.o c_add.o c_sub.o \
cscalar_prod.o cbesi.o cbesip.o cbesjp.o cbesk.o \
cbeskp.o cbesy.o cbesyp.o cbigdet.o cbunk.o ccaldet.o \
cfielder.o cfieldez.o cfieldhr.o cfieldhz.o cfirstcoeff.o cgetrf.o \
cgetri.o cmatmul.o cmatmul2.o cmatsolv.o cmplxtransfer.o cproduct.o \
c_prod.o cscalar_div.o cscalar_div2.o cacon.o cbesh.o cunk1.o cunk2.o \
cgemm.o cgemv.o cswap.o ctrsm.o ctrtri.o cgetf2.o claswp.o \
ctrmm.o ctrti2.o icamax.o cscal.o ctrmv.o cgeru.o

INCLUDE= Makefile f2c.h
# prev, slice, and hist are the executables this file creates
all: tst

# only occasionally necessary for template instantiation problems
clean:
	rm -f *.o *core

# add one rule per executable as follows:
tst: tst.c $(OBJS)
	$(CC) -g -o tst $(OBJS) $(LDFLAGS) $(LIBS)
test: test.c $(OBJS)
	$(CC) -g -o test $(OBJS) $(LDFLAGS) $(LIBS)

# a few C files need extra compiler arguments
tst.o : tst.c $(INCLUDE)
	$(CC) $(CFLAGS) -c tst.c
