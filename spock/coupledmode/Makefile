#
# Waveguide Makefile 2/2/03
#

CC = gcc
LD = ld -r
RM = rm -f
CFLAGS =  -g -c

OBJS= cquanc8.o bgsp.o agsp.o sinus.o

INCLUDE= ./../

libcoupledmode.o : 
	${RM} *.o
	$(CC) -g -c -I${INCLUDE} *.c
	${LD} ${OBJS} -o libcouplemode.o

%.o: %.c
	$(CC) -I${INCLUDE} $(CFLAGS) -c $<
