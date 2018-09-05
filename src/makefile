CC := gcc
CFLAGS :=  -Wall -lm -fopenmp -lpython3.6m -lstdc++ -Wl,--no-as-needed -ldl -lutil
INCLUDES:= -I/usr/include/python3.6m/
#INCLUDES:= -I/home/reemha/anaconda3/include/python3.6m
DEPS := utils.o problem.o initialize.o python_bridge.o
CFILES :=  utils.o problem.o

all: main test

main: main.c $(DEPS) structs.h
	$(CC) $(DEPS) -o main main.c $(CFLAGS)

test: test.c utils.o
	$(CC) utils.o -o test test.c $(CFLAGS)

initialize: structs.h python_bridge.h initialize.c initialize.h
	$(CC) -c -o initialize.o initialize.c $(CFLAGS) $(INCLUDES)

pythonb: python_bridge.h python_bridge.c
	$(CC) -c -o python_bridge.o python_bridge.c $(CFLAGS) $(INCLUDES)

utils: utils.c utils.h
	$(CC) -c -o utils.o utils.c $(CFLAGS)

problem: problem.c problem.h structs.h
	$(CC) -c -o problem.o problem.c $(CFLAGS)

clean:
	rm -f *.o 
