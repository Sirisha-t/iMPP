#ifndef OVERLAPGRAPH_H_
#define OVERLAPGRAPH_H_
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


/*
adjacency_list<         // a 2D structure, the first is a vertex, vertex contains its edge list
   OutEdgeList,         // container is used to represent the edge lists, vecS faster to iterate, listS faster to add
   VertexList,          // container is used to represent the outer two-dimensional container
   Directed,            // directed / undirected / (bidirectional) directed with access to both the in-edges and out-edges
   VertexProperties,    //
   EdgeProperties,      //
   GraphProperties,     //
   EdgeList>            //
*/
class STRVertexType
{
public:
  explicit STRVertexType(const int i)       { rid_ = i;       }
  inline bool IsSeed(void)                  { return ff_seed; }
  /* use for overlap graph */
  int rid_ = -1;                // the read ID for the current vertex
  int len_ = 0;                 // the length of the read
  bool ff_delete = false;       // tag to delete after condense
  bool ff_seed = false;         // tag to stop travering
};

class STREdgeType
{
public:
  inline void SetCondensedTag(const bool i) { ff_condensed = i; }
  inline bool IsCondensed(void)             { return ff_condensed; }
  inline bool IsCycle(void)                 { return ff_cycle; }
  /* use for overlap graph */
  int sid_ = -1;
  bool ff_condensed = false;             // tag to check whether the edge has been condensed
  bool ff_delete = false;                // tag to delete after condense cycle
  bool ff_cycle = false;
  IntegerVector path_info_;
  /* a vector contanins the path information of the edge with the following format:
     read_ID:read_len:p_A: read_ID:read_len:p_A ...
  */
  int p_A = -1;
};

typedef boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, STRVertexType, STREdgeType> BoostSTRGraph;
typedef boost::graph_traits<BoostSTRGraph>::vertex_descriptor BoostSTRVertex;
typedef boost::graph_traits<BoostSTRGraph>::edge_descriptor   BoostSTREdge;
typedef boost::graph_traits<BoostSTRGraph>::vertex_iterator   BoostSTRVertexIter;
typedef boost::graph_traits<BoostSTRGraph>::edge_iterator     BoostSTREdgeIter;

typedef struct Cycle
{
  BoostSTRVertex  head;
  BoostSTRVertex  tail;
  IntegerVector   path;
  bool            ff_cycle;
} cycle_s;

double parseEvalue(const std::string& e);
bool formCycle(IntegerList & path, int new_tail);
void sortAnchor(std::unordered_map<int, anchors>& anchorperf, anchors& anchor_set);
void concateAnchor(std::unordered_map<int, anchors>& anchorperf, anchors& anchor_set);
unsigned int goodMapping(std::string& CIAGR, int flag);
bool goodJunction(BooleanVector& coverage, int pos_i, int pos_j);
bool thisIsOrphat(std::string& header);


class StrGraph {
private:
  bool ff_stringGraph = false;        // false if overlap graph
  bool ff_spadesGraph = false;        // true if the graph from SPAdes
  // the overlap graph
  int             p_node_id = 0;      // both
  int             p_order;
  BoostSTRGraph*  p_graph_;           // both
  StringVector    p_seq;              // both
  IntegerVector   p_all_nodes;
  BooleanVector   p_traversed;
  BooleanVector   p_crossed;
  BooleanVector p_crossed_h;
  BooleanVector p_crossed_t;
  std::unordered_map<int, BoostSTRVertex> vertex; // both

  void markVertexAsTraversed(int v_id);
  void condense(const BoostSTRVertex seed,
    std::list<cycle_s>& to_add_cycle);
  void condense(const BoostSTREdge source_edge,
    std::list<cycle_s>& to_add_cycle);
  std::string RC(int i);

  // for fgs predicted read
  String2Integer readNameMap;
  Integer2Integer predMap;

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

  void DirectedDFS(struct anchor& anchor, predEdges& prededge_set, int MAX_EXTEND_LENGTH, char* path_name, int& pathnum, bool & ff_spadesGraph);
  void savePath(char* path_name, const struct anchor& anchor);
  //void compressPath(char* path_name);

  // SPAdes graph
  String2Integer edgeId;
  int kmer_length;
  void savePath(char* path_name, const IntegerList& path);
  int getPathLength(const IntegerList& path);
public:
  StrGraph();
  ~StrGraph();
  void showGraph();

  // the overlap graph
  bool readFGSPredictions(char* fastaname, char* gffname);
  bool readAsqgFile(char* filename);
  void CondenseGraph();
  void writeGraph(char* filename);

  // the string graph
  bool readSGFile(std::string & filename);
  bool extendGraph(
    char* fgsoutname,
    int MAX_EXTEND_LENGTH,
    char* path_name,
    bool db_flag);
  // SPAdes graph
  bool readSPAdesFile(std::string & filename, int km);
  void extendSPAdes(char* path_name, int MAXLEN = 500);

  //merge graph
  void mergeSG(std::string& sam_name, std::string& contig_name);
  void writeMergedGraph(std::string & filename);

};
#endif
