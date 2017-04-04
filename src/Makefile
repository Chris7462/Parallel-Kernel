CPP=g++
CPPFLAGS=-O3
OMPFLAG=-fopenmp
PGFLAG=-pg

default:
	$(CPP) $(CPPFLAGS) main.cc -o main

omp:
	$(CPP) $(CPPFLAGS) $(OMPFLAG) main.cc -o main

pg:
	$(CPP) $(CPPFLAGS) $(PGFLAG) $(OMPFLAG) main.cc -o main

clean:
	rm main *.out
