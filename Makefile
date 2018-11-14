MF=	Makefile

CC=	mpicc
CFLAGS= -O3 -g -Wall

LFLAGS=	-lm 

EXE=	build/image

SRC= \
	src/imagenew.c \
	src/pgmio.c\
	src/arralloc.c

INC=\
	include/pgmio.h\
	include/arralloc.h
OBJ=\
	build/imagenew.o\
	build/pgmio.o\
	build/arralloc.o

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .o

$(EXE):$(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

#OBJ=	$(SRC:.c=.o)
build/%.o: src/%.c 
	$(CC) $(CFLAGS) -c $< -o $@
	 
#.c.o:
#	$(CC) $(CFLAGS) -c $<

#all:	$(EXE)

#$(EXE):	$(OBJ)
#	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

#$(OBJ):	$(INC)

#$(OBJ):	$(MF)

clean:
	rm -f $(OBJ) $(EXE) 
