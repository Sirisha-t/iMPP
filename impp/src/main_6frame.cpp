#include "six_frame.h"

int main(int argc, char* argv[])
{ 
  int c;
  char* fgspred = (char*)malloc(sizeof(char)*100);
  char* bwaedge = (char*)malloc(sizeof(char)*100);
  char* bwapath = (char*)malloc(sizeof(char)*100);
  char* fastaread = (char*)malloc(sizeof(char)*100);
  char* outf = (char*)malloc(sizeof(char)*100); 
 //std::string graph_name = "", path_out = "", gff_name = "";
  while ((c = getopt(argc, argv, "f:e:p:r:o:")) >= 0)
  {
    if (c == 'f')	fgspred = optarg;
    else if(c == 'e')	bwaedge = optarg;
    else if(c == 'p') 	bwapath = optarg;
    else if(c == 'r') 	fastaread = optarg;
    else if(c == 'o') 	outf = optarg;

    else           std::cout<<"ERROR: Check input arguments! \n";
  }

  if(fgspred == "" || bwaedge == "" || bwapath == ""){
    std::cout<<"ERROR: Missing input arguments. See usage options below \n";
  }
  else{
  FGC* p_s = new FGC(); 
   p_s->loadRead(fgspred);
   p_s->loadBWA(bwaedge);
   p_s->loadBWA(bwapath);
   p_s->OutSeq(fastaread, outf);

  
  //delete p_g;
  }
  //free(graph_name);
  return 0;
}

