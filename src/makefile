all: soNumerit soEF so

main: main.c 2d_kershaw.o utils.o
	gcc -o main main.c -lm -fopenmp

soEF: SuOlsonEddingtonFactor.c
	gcc -o soEF SuOlsonEddingtonFactor.c -lm -fopenmp

so: tridFunc.o SuOlsonAll.c SuOlsonAll.h
	gcc -o so SuOlsonAll.c tridFunc.o -lm -fopenmp

utils: utils.c utils.h
	gcc -c -o utils.o utils.c -lm -fopenmp

clean:
	rm -f *.o soNumerit
