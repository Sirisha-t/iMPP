CC=g++
CFLAG= -std=c++11 -Wall -Wno-deprecated -fopenmp
OPT = -O3

all: SG MERGE IMPP ASQG

SG: string_graph.o main_sg.o
	$(CC) $(CFLAG) $(OPT) -o ../bin/impp_getSG main_sg.o string_graph.o 
MERGE: merge_sg.o extend_anchor.o unionFind.o
	$(CC) $(CFLAG) $(OPT) -o ../bin/impp_mergeSG extend_anchor.o merge_sg.o unionFind.o 
IMPP: main_anchor.o extend_anchor.o unionFind.o
	 $(CC) $(CFLAG) $(OPT) -o ../bin/impp_extendAnchor main_anchor.o extend_anchor.o unionFind.o
ASQG: main_AssemblyGraph.o AssemblyGraph.o GraphEssential.o GraphPrune.o StringUtils.o 
	$(CC) $(CFLAG) $(OPT) -o ../bin/impp_modASQG main_AssemblyGraph.o AssemblyGraph.o GraphEssential.o GraphPrune.o StringUtils.o
main_AssemblyGraph.o: main_AssemblyGraph.cc AssemblyGraph.h 
	$(CC) $(CFLAG) $(OPT) -c -o main_AssemblyGraph.o main_AssemblyGraph.cc
AssemblyGraph.o: AssemblyGraph.cc AssemblyGraph.h
	$(CC) $(CFLAG) $(OPT) -c -o AssemblyGraph.o AssemblyGraph.cc 
GraphEssential.o: GraphEssential.cc GraphEssential.h GraphNodeType.h GraphEdgeType.h
	$(CC) $(CFLAG) $(OPT) -c -o GraphEssential.o GraphEssential.cc
GraphPrune.o: GraphPrune.cc GraphPrune.h DataType.h
	$(CC) $(CFLAG) $(OPT) -c -o GraphPrune.o GraphPrune.cc
StringUtils.o: StringUtils.cc StringUtils.h 
	$(CC) $(CFLAG) $(OPT) -c -o StringUtils.o StringUtils.cc
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
