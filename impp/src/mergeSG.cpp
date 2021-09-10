#include <unistd.h>
#include "SG.h"
int mergeSG_usage()
{
  std::cout << "\n";
  std::cout << "Usage:   impp mergeSG <string_graph.fq> <bwa.sg.contigs.sam> <contigs.fq> \n\n";
  std::cout << "         string_graph.fq      string graph after calling assembly\n";
  std::cout << "         bwa.sg.contigs.sam   bwa output with the mapping of <string_graph.fq> to <contigs.fq>\n";
  std::cout << "         contigs.fq           contig file used to connect edges in <string_graph.fq>\n";
  std::cout << "\n";
	return 1;
}

int main(int argc, char* argv[])
{
  if(argc != 4) return mergeSG_usage();
  std::string sg_name = argv[1];
  std::string sam_name = argv[2];
  std::string contig_name = argv[3];
  std::cout << "Loading string graph...\n";
  StrGraph* p_g = new StrGraph();
  p_g->readSGFile(sg_name);
  std::cout << "Merging string graph...\n";
  p_g->mergeSG(sam_name, contig_name);
  std::cout << "Writing string graph...\n";
  p_g->writeMergedGraph(sg_name);
  delete p_g;
  return 0;
}
