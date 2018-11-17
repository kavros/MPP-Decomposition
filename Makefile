MF=	Makefile

CC=	mpicc
CFLAGS= -O3 -g -Wall

LFLAGS=	-lm 

EXE=	build/image

SRC= \
	src/imagenew_v2.c \
	src/pgmio.c\
	src/arralloc.c

INC=\
	include/pgmio.h\
	include/arralloc.h
OBJ=\
	build/imagenew_v2.o\
	build/pgmio.o\
	build/arralloc.o

nThreads=1
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
	mpirun -n ${nThreads} ./build/image

qsub:
	qsub scripts/image.pbs

validation:
	diff data/output/imagenew192x128_expectedOutput.pgm  data/output/imagenew192x128.pgm
