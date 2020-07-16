#include "SG.h"

int main(int argc, char* argv[])
{ 
  int c, max_length = 150;
  char* graph_name = (char*)malloc(sizeof(char)*100);
  char* path_out = (char*)malloc(sizeof(char)*100);
  char* gff_name = (char*)malloc(sizeof(char)*100);
  //std::string graph_name = "", path_out = "", gff_name = "";
  while ((c = getopt(argc, argv, "g:a:l:p:")) >= 0)
  {
    if (c == 'g')  graph_name = optarg;
    else if (c == 'a')  gff_name = optarg;
    else if (c == 'l')  max_length = std::stoi(optarg);
    else if (c == 'p')  path_out = optarg;
    else                std::cout<<"ERROR: Check input arguments! \n";
  }

  if(graph_name == ""){
    std::cout<<"ERROR: Please provide input graph file (option -g). See usage options below \n";
  }
  else{
  StrGraph* e_g = new StrGraph();
  e_g->readFqFile(graph_name);
  e_g->DirectedDFS(gff_name, max_length, path_out );
  
  delete e_g;
  }
  //free(graph_name);
  //free(path_out);
  //free(gff_name);

  return 0;
}
