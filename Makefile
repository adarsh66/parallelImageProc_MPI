
MF=	Makefile

#Select DNDEBUG if you do not want to see debug msgs
debug= -DNDEBUG
#debug=

#morar
CC=	mpicc
DEBUG ?= 1
CFLAGS = -fastsse $(debug)

#personal laptop
#CC=	cc
#CFLAGS= -g

LFLAGS=	-lm

EXE	=	imagerecreation

SRC= \
	imagerecreation.c \
	src/inputfuncs.c \
	src/pgmio.c 

INC=\
	src/pgmio.h \
	src/inputfuncs.h \
	src/log.h
#
# No need to edit below this line
#


.SUFFIXES:
.SUFFIXES: .c .o

OBJ=	$(SRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -o $(<:.c=.o) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(INC)

$(OBJ):	$(MF)

clean:
	rm -f $(OBJ) $(EXE) core
