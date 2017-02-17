CC=gcc
CFLAGS=-Wall -pipe -pedantic
INC =   -I /usr/include
LDFLAGS = -lm 

OBJS= cdiv.o cexp.o csqrt.o \
      clog.o algaas2_matsys.o

INCLUDE= Makefile util.h
# algaas2 is the executable this makefile creates
all: algaas2

# only occasionally necessary
clean:
	rm -f *.o *core

# add one rule per executable as follows:
algaas2: algaas2_matsys.c $(OBJS)
	$(CC) -o algaas2 $(OBJS) $(LDFLAGS)

