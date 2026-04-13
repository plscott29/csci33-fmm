all: fmm direct compare data_writer

fmm:
	gcc-15 -Wall -Wno-unused-function -std=c99 -fopenmp -o fmm fmm.c -lm

direct:
	gcc-15 -Wall -Wno-unused-function -std=c99 -fopenmp -o direct direct.c -lm

compare:
	gcc-15 -Wall -Wno-unused-function -std=c99 -o compare compare.c -lm

data_writer:
	gcc-15 -Wall -Wno-unused-function -std=c99 -o data_writer data_writer.c -lm

clean:
	rm -f fmm direct compare data_writer