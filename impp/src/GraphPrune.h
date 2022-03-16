//  Author: Cuncong Zhong
//  Last modification: 11/10/2021

#ifndef __GRAPHPRUNE_H_
#define __GRAPHPRUNE_H_

#include <queue>

#include "StringUtils.h"
#include "GraphEssential.h"

class GraphPrune    {
  public:
    // void construction function, does nothing
    GraphPrune()    {}
    // void destruction function, does nothing
    ~GraphPrune()   {}

    // Remove orphan vertices
    // parameter:
    //    g: the pointer to the graph where the function operates on
    void RemoveOrphanVertices(AssemblyGraphType *g);

    // Condence single paths
    // parameter:
    //    g: the pointer to the graph where the function operates on
    void CondenceSinglePaths(AssemblyGraphType *g);

    // Resolve read orientation; for each read only one direction is retained in the assembly graph
    // Algorithm: Start with an unresolved node with the highest degree and set it to positive strand
    // Repeat:  set its neigbhors' orientation according to the edge;
    //          if the neighbor has already been set and is conflicting with the current edge, remove the current edge
    //          mark the neighbors as resolved
    //          return to step 1
    // Internal data structure used: priority_queue
    // parameter:
    //    g: the pointer to the graph where the function operates on
    void ResolveOrientation(AssemblyGraphType *g);

    // Resolve the sequences in the assembly graph.
    // Use reverse complementary if the read is assigned negative strand
    // Also update all edges as no reverse-completementary
    // Parameter:
    //    g: the pointer to the graph where the function operates on
    void ResoveSequence(AssemblyGraphType *g);  
  
  private:
    
};

#endif          // __GRAPHPRUNE_H_
