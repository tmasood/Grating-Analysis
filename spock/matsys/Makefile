#
# Waveguide Makefile 2/2/03
#

CC = gcc
LD = ld -r
RM = rm -f

OBJS= ingaasp_adachi.o algainas_adachi.o algainas_mondry.o ingaasp_henry.o

INCLUDE= ./../

libmatsys.o : 
	${RM} *.o
	$(CC) -g -c -I${INCLUDE} *.c
	${LD} ${OBJS} -o libmatsys.o
