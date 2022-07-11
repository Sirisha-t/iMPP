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
#include <math.h>  

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


class StrGraph {
private:
  bool ff_stringGraph = false;        // false if overlap graph
  bool ff_spadesGraph = false;        // true if the graph from SPAdes
  // the overlap graph
 // int             p_node_id = 0;      // both
  int             p_order;
 // BoostSTRGraph*  p_graph_;           // both
 // StringVector    p_seq;              // both
  IntegerVector   p_all_nodes;
  BooleanVector   p_traversed;
  BooleanVector   p_crossed;
 // std::unordered_map<int, BoostSTRVertex> vertex; // both

  void markVertexAsTraversed(int v_id);
  void condense(const BoostSTRVertex seed,
    std::list<cycle_s>& to_add_cycle);
  void condense(const BoostSTREdge source_edge,
    std::list<cycle_s>& to_add_cycle);
  std::string RC(int i);

public:
  StrGraph();
  ~StrGraph();
  void showGraph();

  int             p_node_id = 0;      // both
  int             p_order;
  BoostSTRGraph*  p_graph_;           // both
  StringVector    p_seq;              // both
  std::unordered_map<int, BoostSTRVertex> vertex; // both

  // the overlap graph
  bool readAsqgFile(char* filename);
  void CondenseGraph();
  void writeGraph(char* filename);
};
#endif
