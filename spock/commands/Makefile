#
# Waveguide Makefile 2/2/03
#

CC = gcc
LD = ld -r
RM = rm -f

OBJS= createlayerfile.o genindexprofile.o genplotindexprofile.o

INCLUDE= ./../

libcommands.o : 
	${RM} *.o
	$(CC) -g -c -I${INCLUDE} *.c
	${LD} ${OBJS} -o libcommands.o
