CC := gcc
CFLAGS :=  -Wall -lm -fopenmp
DEPS := utils.o problem.o initialize.o
CFILES :=  utils.o problem.o

all: main

main: main.c $(DEPS) structs.h
	$(CC) $(DEPS) -o main main.c $(CFLAGS)

initialize: structs.h initialize.c initialize.h
	$(CC) -c -o initialize.o initialize.c $(CFLAGS)

utils: utils.c utils.h
	$(CC) -c -o utils.o utils.c $(CFLAGS)

problem: problem.c problem.h structs.h
	$(CC) -c -o problem.o problem.c $(CFLAGS)

clean:
	rm -f *.o 
