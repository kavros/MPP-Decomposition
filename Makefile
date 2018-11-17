MF=	Makefile

CC=	mpicc
CFLAGS= -O3 -g -Wall

LFLAGS=	-lm 

EXE=	build/image

SRC= \
	src/imagenew_v2.c \
	src/pgmio.c\
	src/arralloc.c\
	src/argtable3.c

INC=\
	include/pgmio.h\
	include/arralloc.h\
	include/argtable3.h
OBJ=\
	build/imagenew_v2.o\
	build/pgmio.o\
	build/arralloc.o\
	build/argtable3.o

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

qsub:
	qsub scripts/image.pbs

validation:
	diff data/output/imagenew192x128_expectedOutput.pgm  data/output/imagenew192x128.pgm
