CC := gcc
CFLAGS := -fopenmp -I/usr/include/python2.7 -I/usr/include/x86_64-linux-gnu/python2.7  -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security  -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -I.
LDFLAGS := -lm -L/usr/lib/python2.7/config-x86_64-linux-gnu -L/usr/lib -lpython2.7 -lpthread -ldl  -lutil -lm  -Xlinker -export-dynamic -Wl,-O1 -Wl,-Bsymbolic-functions
DEPS := utils.o problem.o initialize.o

all: main

main: main.c $(DEPS) structs.h
	$(CC) $(CFLAGS) $(DEPS) main.c -o main.o $(LDFLAGS)

initialize: structs.h initialize.c initialize.h
	$(CC) $(CFLAGS) initialize.c -o initialize.o $(LDFLAGS)  

utils: utils.c utils.h
	$(CC) -c -o utils.o utils.c $(CFLAGS)

problem: problem.c problem.h structs.h
	$(CC) -o problem.o problem.c $(LDFLAGS) $(CFLAGS) 

clean:
	rm -f *.o 
