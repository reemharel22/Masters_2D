CC := @gcc
CFLAGS := -fopenmp -g -O2 -w
LDFLAGS := -lm -O2
DEPS := utils.o problem.o initialize.o calculations.o diagnostics.o

all: utils initialize problem calculation diagnostics main

main: main.c $(DEPS) structs.h
	$(CC) $(CFLAGS) $(DEPS) main.c -o main.o $(LDFLAGS)

initialize: structs.h initialize.c initialize.h
	$(CC) $(CFLAGS) initialize.c -c -o initialize.o $(LDFLAGS)  

calculation: calculations.c calculations.h structs.h
	$(CC) $(CFLAGS) calculations.c -c -o calculations.o -lm -O2

utils: utils.c utils.h
	$(CC) $(CFLAGS) utils.c -c -o utils.o $(CFLAGS)

diagnostics: diagnostics.c diagnostics.h
	$(CC) $(CFLAGS) -c diagnostics.c -o diagnostics.o -O2

problem: problem.c problem.h structs.h
	$(CC) $(CFLAGS) problem.c -c -o problem.o $(LDFLAGS) 

test: test.c utils.o
	$(CC) $(CFLAGS) utils.o test.c -o test.o

clean:
	rm -f *.o 
