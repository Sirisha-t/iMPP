CC=g++
CFLAG= -std=c++11 -Wall -Wno-deprecated -fopenmp
OPT = -O3

all: EVAL PEPSTATS
EVAL: eval.o
	$(CC) $(CFLAG) $(OPT) -o ../bin/eval eval.o  
PEPSTATS: pepstats.o
	$(CC) $(CFLAG) $(OPT) -o ../bin/pepstats pepstats.o
eval.o: eval.cpp eval.h
	$(CC) $(CFLAG) $(OPT) -c -o eval.o eval.cpp
pepstats.o:
	$(CC) $(CFLAG) $(OPT) -c -o pepstats.o peptideStats.cpp 

clean:
	rm -rf *.o ../bin/* *~

