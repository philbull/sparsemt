
all: sparse.o
	gcc ghmock.o -lm -lgsl -lgslcblas -lgomp -o ghmock

sparse.o: sparse.c
	gcc -O4 -fopenmp -std=c99 -DSTANDALONE -c ghmock.c

python: 
	python setup.py build_ext --inplace

clean:
	rm -f *.o
	rm -f sparse.so
	rm -rf build/
