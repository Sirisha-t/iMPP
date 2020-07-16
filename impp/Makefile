CC=g++ -c -o
CFLAG= -O3
LINK= g++ -o

all: IMPP

IMPP: 	main.o SG.o
	$(LINK) ../bin/run_impp main.o SG.o -fopenmp
SG.o: 	SG.cpp
	$(CC) SG.o SG.cpp -fopenmp
main.o: main.cpp
	$(CC) main.o main.cpp -fopenmp

clean:
	rm -rf *.o ../bin/* *~
