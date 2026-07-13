all: NpT2DMall

CC=gcc
OPT2=-Wall -Wextra -std=c99 -O3 -march=native -flto -ffast-math -fno-finite-math-only -fno-unsafe-math-optimizations -fno-align-functions -fno-align-jumps -fno-align-loops -fno-align-labels
OPT1=-Wall -std=c99 -g
LIBS3= -lgsl -lgslcblas -lm

NpT2DMall: NpT_2D_M.c
	${CC} ${OPT2} ${LIBS3} -c lists.c
	${CC} ${OPT2} ${LIBS3} -c NpT_2D_M.c
	${CC} ${OPT2} ${LIBS3} -c distance.c
	${CC} ${OPT2} ${LIBS3} -c init.c
	${CC} ${OPT2} ${LIBS3} -c mc_steps.c
	${CC} ${OPT2} ${LIBS3} -c printdata.c
	${CC} ${OPT2} ${LIBS3} -c verletlist.c
	${CC} ${OPT2} ${LIBS3} -c inputdata.c
	${CC} ${OPT2} ${LIBS3} -c printdisplay.c
	${CC} ${OPT2} ${LIBS3} lists.o distance.o init.o mc_steps.o printdata.o verletlist.o inputdata.o printdisplay.o NpT_2D_M.o -o NpT_2D_M.out ${LIBS3}
	
clean:
	rm -f lists.o
	rm -f NpT_2D_M.o
	rm -f distance.o
	rm -f init.o
	rm -f mc_steps.o
	rm -f printdata.o
	rm -f verletlist.o
	rm -f inputdata.o
	rm -f printdisplay.o
	rm -f NpT_2D_M.out
