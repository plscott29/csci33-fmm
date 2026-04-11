all: fmm direct

fmm:
	gcc-15 -Wall -std=c99 -fopenmp -o fmm fmm.c -lm

direct:
	gcc-15 -Wall -std=c99 -fopenmp -o direct direct.c -lm

clean:
	rm -f fmm direct