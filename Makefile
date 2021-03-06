MF=	Makefile

CC=	mpicc
CFLAGS= -O3 -g -Wall

LFLAGS=	-lm 

EXE=	build/image

SRC= \
	src/imagenew_v2.c \
	src/pgmio.c\
	src/arralloc.c\
	src/argtable3.c\
	src/functions.c\
	src/initializations.c\
	src/communications.c

INC=\
	include/pgmio.h\
	include/arralloc.h\
	include/argtable3.h\
	include/functions.h\
	include/initialization.h\
	include/globalVariables.h\
	include/communications.h
	
OBJ=\
	build/imagenew_v2.o\
	build/pgmio.o\
	build/arralloc.o\
	build/argtable3.o\
	build/functions.o\
	build/initializations.o\
	build/communications.o

n=1
#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .o

$(EXE):$(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

build/%.o: src/%.c 
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXE) 

run: 
	mpirun -n ${n} ./build/image

run_delta:
	mpirun -n ${n} ./build/image -d 1

qsub:
	qsub -v MPISIZE=${n},INPUT_IMG="./data/input/edgenew192x128.pgm",OUTPUT_IMG="./data/qsub/imagenew192x128.pgm"   scripts/image.pbs

validation:
	diff data/output/imagenew192x128_expectedOutput.pgm  data/output/imagenew192x128.pgm
