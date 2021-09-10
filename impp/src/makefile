CC=g++
CFLAG= -std=c++11 -Wall -Wno-deprecated -fopenmp
OPT = -O3

all: SG MERGE IMPP

SG: string_graph.o main_sg.o
	$(CC) $(CFLAG) $(OPT) -o ../bin/impp_getSG main_sg.o string_graph.o 
MERGE: merge_sg.o extend_anchor.o unionFind.o
	$(CC) $(CFLAG) $(OPT) -o ../bin/impp_mergeSG extend_anchor.o merge_sg.o unionFind.o 
IMPP: main_anchor.o extend_anchor.o unionFind.o
	 $(CC) $(CFLAG) $(OPT) -o ../bin/impp_extendAnchor main_anchor.o extend_anchor.o unionFind.o
merge_sg.o: mergeSG.cpp SG.h unionFind.hpp
	$(CC) $(CFLAG) $(OPT) -c -o merge_sg.o mergeSG.cpp
string_graph.o: string_graph.cpp SG.h
	$(CC) $(CFLAG) $(OPT) -c -o string_graph.o string_graph.cpp 
main_sg.o: main_SG.cpp SG.h
	$(CC) $(CFLAG) $(OPT) -c -o main_sg.o main_SG.cpp
main_anchor.o: main_extendAnchor.cpp SG.h
	$(CC) $(CFLAG) $(OPT) -c -o main_anchor.o main_extendAnchor.cpp
extend_anchor.o: extend_anchor.cpp SG.h unionFind.hpp
	$(CC) $(CFLAG) $(OPT) -c -o extend_anchor.o extend_anchor.cpp
unionFind.o: unionFind.cpp unionFind.hpp
	$(CC) $(CFLAG) $(OPT) -c -o unionFind.o unionFind.cpp


clean:
	rm -rf *.o ../bin/impp_* *~
