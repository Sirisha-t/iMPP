#include "SG.h"

int main(int argc, char* argv[])
{
  int c, max_length = 150;
  int kmer_len = 0;
  bool ff_graph = false;
  std::string graph_name;
  //char* graph_name = (char*)malloc(sizeof(char)*100);
  char* path_out = (char*)malloc(sizeof(char)*100);
  char* gff_name = (char*)malloc(sizeof(char)*100);
  //std::string graph_name = "", path_out = "", gff_name = "";
  while ((c = getopt(argc, argv, "g:a:l:p:k:")) >= 0)
  {
    if (c == 'g')  graph_name = optarg;
    else if (c == 'a')  gff_name = optarg;
    else if (c == 'l')  max_length = std::stoi(optarg);
    else if (c == 'p')  path_out = optarg;
    else if (c == 'k')  kmer_len = std::stoi(optarg);
    else                std::cout<<"ERROR: Check input arguments! \n";
  }

  if(graph_name == ""){
    std::cout<<"ERROR: Please provide input graph file (option -g). See usage options below \n";
  }
  else{
      StrGraph* e_g = new StrGraph();
      //if(kmer_len > 0){
       // e_g->readSPAdesFile(graph_name, kmer_len);
       // ff_graph = true;
      //}
     // else{
      e_g->readSGFile(graph_name);

      e_g->extendGraph(gff_name, max_length, path_out, ff_graph);

      delete e_g;
    }

  return 0;
}
