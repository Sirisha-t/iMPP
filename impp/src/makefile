CC=g++
CFLAG= -std=c++11 -Wall -Wno-deprecated -fopenmp
OPT = -O3

all: SG IMPP

SG: string_graph.o main_sg.o
	$(CC) $(CFLAG) $(OPT) -o ../bin/impp_getSG main_sg.o string_graph.o 
IMPP: main_anchor.o extend_anchor.o
	 $(CC) $(CFLAG) $(OPT) -o ../bin/impp_extendAnchor main_anchor.o extend_anchor.o 
#FGC_6f: main_6f.o six_frame.o
#	$(CC) $(CFLAG) $(OPT) -o ../bin/FGC_get6f main_6f.o six_frame.o 
string_graph.o: string_graph.cpp SG.h
	$(CC) $(CFLAG) $(OPT) -c -o string_graph.o string_graph.cpp 
main_sg.o: main_SG.cpp SG.h
	$(CC) $(CFLAG) $(OPT) -c -o main_sg.o main_SG.cpp
main_anchor.o: main_extendAnchor.cpp SG.h
	$(CC) $(CFLAG) $(OPT) -c -o main_anchor.o main_extendAnchor.cpp
extend_anchor.o: extend_anchor.cpp SG.h
	$(CC) $(CFLAG) $(OPT) -c -o extend_anchor.o extend_anchor.cpp
#main_6f.o: main_6frame.cpp six_frame.h
#	$(CC) $(CFLAG) $(OPT) -c -o main_6f.o main_6frame.cpp
#six_frame.o: six_frame.cpp six_frame.h
#	$(CC) $(CFLAG) $(OPT) -c -o six_frame.o six_frame.cpp


clean:
	rm -rf *.o ../bin/impp_* *~
