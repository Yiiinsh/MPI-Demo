CC = mpicc
CFLAGS = -std=c99 -O3

EXE = solution 

SRC = main.c pgmio.c arraytool.c
INC = pgmio.h arraytool.h
OBJ = $(SRC:.c=.o)

all : $(EXE)

$(EXE) : $(OBJ)
	   $(CC) $(CFLAGS) -o $@ $^

$(OBJ) : %o : %c $(INC)
	$(CC) $(CFLAGS) -c $<


.PHONY: clean
clean :
	rm -f $(OBJ) $(EXE)

