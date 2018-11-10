MF=	Makefile

CC=	mpicc
CFLAGS= -O3 -g

LFLAGS=	-lm

EXE=	image

SRC= \
	imagenew.c \
	pgmio.c

INC=\
	pgmio.h

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .o

OBJ=	$(SRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(INC)

$(OBJ):	$(MF)

clean:
	rm -f $(OBJ) $(EXE) core
