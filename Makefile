CC = mpicc
CFLAGS = -std=c99 -O3

PARALLEL_EXE = parallel 
SERIAL_EXE = serial

SRC = main.c pgmio.c arraytool.c boundary.c
INC = pgmio.h arraytool.h defs.h imgprocessing.h boundary.h
PARALLEL_SRC = parallel_imgprocessing.c
SERIAL_SRC = serial_imgprocessing.c
OBJ = $(SRC:.c=.o)
PARALLEL_OBJ = parallel_imgprocessing.o
SERIAL_OBJ = serial_imgprocessing.o

all : $(PARALLEL_EXE)

$(PARALLEL_EXE) : $(OBJ) $(PARALLEL_OBJ)
	$(CC) $(CFLAGS) -o $@ $^

$(OBJ) : %o : %c $(INC)
	$(CC) $(CFLAGS) -c $<

.PHONY: serial
serial : $(SERIAL_EXE)

$(SERIAL_EXE) : $(OBJ) $(SERIAL_OBJ)
	$(CC) $(CFLAGS) -o $@ $^

$(PARALLEL_OBJ) : $(PARALLEL_SRC) $(INC)
	$(CC) $(CFLAGS) -c $<

$(SERIAL_OBJ) : $(SERIAL_SRC) $(INC)
	$(CC) $(CFLAGS) -c $<

.PHONY: test
test :
	./test
	
.PHONY: clean
clean :
	rm -f $(OBJ) $(PARALLEL_EXE) $(SERIAL_EXE) $(PARALLEL_OBJ) $(SERIAL_OBJ)

