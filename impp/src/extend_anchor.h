#ifndef OVERLAPGRAPH_H_
#define OVERLAPGRAPH_H_
#include "string_graph.h"
#include <boost/tokenizer.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <list>
#include <stack>
#include <time.h>
#include <math.h>       // log10

typedef std::list<int>                        IntegerList;
struct anchor
{
  IntegerList path;
  int M_upstream = -1;
  int M_downstream = -1;
  /*int RF = -1;
  double E = -1; */
};
struct pEdgeInfo
{
  int node_id;
  int start = -1;
  int end = -1;
  int dir =-1;
  int frame = -1;

};
/*struct edgeP //new with only one struct
{
  IntegerList start;
  IntegerList end;
  IntegerList dir;
  IntegerList frame;

};*/

bool compare_nocase (const struct anchor& first, const struct anchor& second);

typedef std::unordered_map<int, struct pEdgeInfo>  predEdges;
typedef std::list<struct anchor>              anchors;
typedef std::vector<int>                      IntegerVector;
typedef std::vector<IntegerVector>            IntegerVector2D;
typedef std::vector<bool>                     BooleanVector;
typedef std::vector<std::string>              StringVector;
typedef std::unordered_map<int, std::string>  Integer2String;
typedef std::unordered_map<std::string, int>  String2Integer;
typedef std::unordered_map<int, int>          Integer2Integer;
typedef std::set<int>                         IntegerSet;


double parseEvalue(const std::string& e);
bool formCycle(IntegerList & path, int new_tail);
void sortAnchor(std::unordered_map<int, anchors>& anchorperf, anchors& anchor_set);
void concateAnchor(std::unordered_map<int, anchors>& anchorperf, anchors& anchor_set);

class ExtendGraph : public StrGraph {
private:

  // the string graph
  Integer2Integer p_read_length;
  IntegerVector2D p_read_on_node;
  IntegerVector2D p_pos_on_node;
  bool loadAnchor(
    char* uEdgeSet,
    int MAX_EXTEND_LENGTH,
    anchors& anchor_set);
  bool parseAnchorLine(
    int & node,
    struct anchor& new_anchor,
    int MAX_EXTEND_LENGTH);

  // FGS predicted predicted edges
    Integer2Integer pEdge;
/* Functions to load the predicted edges from fgs */
  bool loadPredictedEdge(
    char* gffname,
    predEdges& pEdge_set);
  int parseGFFLine(
    const std::string& line,
    struct pEdgeInfo& new_edge);

  void DirectedDFS(struct anchor& anchor, predEdges& prededge_set, int MAX_EXTEND_LENGTH, char* path_name, int& pathnum, bool &ff_spadesGraph);
  void savePath(char* path_name, const struct anchor& anchor);
  //void compressPath(char* path_name);

  // SPAdes graph
  String2Integer edgeId;
  int kmer_length;
  void savePath(char* path_name, const IntegerList& path);
  int getPathLength(const IntegerList& path);
public:
   ExtendGraph();
  ~ExtendGraph();

  // the string graph
  bool readSGFile(char* filename);
  bool extendGraph(
    char* fgsoutname,
    int MAX_EXTEND_LENGTH,
    char* path_name,
    bool  db_flag);
  // SPAdes graph
  bool readSPAdesFile(char* filename, int km);
  void extendSPAdes(char* path_name, int MAXLEN = 300);
};
#endif
