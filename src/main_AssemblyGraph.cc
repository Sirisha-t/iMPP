//  Author: Cuncong Zhong
//  Last modification: 11/09/2021

#include <iostream>

#include "AssemblyGraph.h"

using namespace std;

int main(int argc, char* argv[])  {
    
    std::string asqg_graph;
    int c;
    clock_t t;
    //char* asqg_graph = (char*)malloc(sizeof(char)*100);
    if((c = getopt(argc, argv, "g:")) >=0){
	if(c == 'g'){
		asqg_graph = optarg;
	}
	else std::cout<<"ERROR: Please check input arguments\n";
    }

    AssemblyGraph *graph = new AssemblyGraph;
    if( graph->IsInitialized()){
    string outname = asqg_graph.substr(0,asqg_graph.find_last_of('.')) + ".mod.asqg";
    cout<<"Graph initialized\n";
    t=clock();
    graph -> RenameASQG(asqg_graph);
    cout<<"DEBUG: after renaming asqg"<<endl;
    string asqg_rename = asqg_graph.substr(0,asqg_graph.find_last_of('.')) + ".rename.asqg";
    graph->LoadGraphASQG(asqg_rename);   
    cout << "DEBUG: after loading graph:" << endl;
    graph->CheckGraphValidity();
    graph->PrintInfo(false);
    t= clock() - t;
    cout<<"Time to load asqg : "<<((double)t)/CLOCKS_PER_SEC<<" s\n";     

    t=clock();
    graph->RemoveOrphanVertices();
    cout << "DEBUG: after removing orphan reads:" << endl;
    graph->CheckGraphValidity();
    graph->PrintInfo(false);
    t= clock() - t;
    cout<<"Time to remove orphan reads : "<<((double)t)/CLOCKS_PER_SEC<<" s\n"; 
    
    t=clock();
    graph->ResolveOrientation();
    cout << "DEBUG: after resolving orientation:" << endl;
    graph->CheckGraphValidity();
    graph->PrintInfo(false);
    t= clock() - t;
    cout<<"Time to resolve orientation : "<<((double)t)/CLOCKS_PER_SEC<<" s\n"; 

    t=clock();
    graph->ResolveSequence();
    cout << "DEBUG: after resolving sequence:" << endl;
    graph->CheckGraphValidity();
    graph->PrintInfo(false);
    t= clock() - t;
    cout<<"Time to resolve sequence : "<<((double)t)/CLOCKS_PER_SEC<<" s\n"; 
    

    t=clock();
    graph->WriteGraphASQG(outname);
    t= clock() - t;
    cout<<"Time to write asqg: "<<((double)t)/CLOCKS_PER_SEC<<" s\n"; 

}

    return 0;
}
