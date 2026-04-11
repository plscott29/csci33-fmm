all: fmm direct compare

fmm:
	gcc-15 -Wall -std=c99 -fopenmp -o fmm fmm.c -lm

direct:
	gcc-15 -Wall -std=c99 -fopenmp -o direct direct.c -lm

compare:
	gcc-15 -Wall -Wno-unused-function -std=c99 -o compare compare.c -lm

clean:
	rm -f fmm direct compare