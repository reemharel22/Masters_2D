CC := gcc
CFLAGS := -fopenmp -I/usr/include/python2.7 -I/usr/include/x86_64-linux-gnu/python2.7  -fno-strict-aliasing -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security  -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -I.
LDFLAGS := -lm -L/usr/lib/python2.7/config-x86_64-linux-gnu -L/usr/lib -lpthread -ldl  -lutil -lm  -Xlinker -export-dynamic -Wl,-O1 -Wl,-Bsymbolic-functions
DEPS := utils.o problem.o initialize.o calculations.o diagnostics.o

all: main

main: main.c $(DEPS) structs.h
	$(CC) $(CFLAGS) $(DEPS) main.c -o main.o $(LDFLAGS)

initialize: structs.h initialize.c initialize.h
	$(CC) $(CFLAGS) initialize.c -o initialize.o $(LDFLAGS)  

calculation: calculations.c calculations.h struct.h
	$(CC) calculations.c -o calcualtions.o -lm -O2

utils: utils.c utils.h
	$(CC) -c -o utils.o utils.c $(CFLAGS)

diagnostics: diagnostics.c diagnostics.h
	$(CC) diagnostics.c -o diagnostics.o -O2

problem: problem.c problem.h structs.h
	$(CC) -o problem.o problem.c $(LDFLAGS) $(CFLAGS) 

test: test.c utils.o
	$(CC) $(CFLAGS) utils.o test.c -o test.o

clean:
	rm -f *.o 
