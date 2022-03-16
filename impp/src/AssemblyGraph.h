//  Author: Cuncong Zhong
//  Last modification: 11/10/2021

#ifndef __ASSEMBLYGRAPH_H_
#define __ASSEMBLYGRAPH_H_

#include "GraphEssential.h"
#include "GraphPrune.h"

class AssemblyGraph : public GraphEssential {
  public:
    // constructor (defined in GraphEssential): initialize *graph_ptr_ 
    // destructor (defined in GraphEssentual): if initialized, delete *graph_ptr

    // wrapper function to call GraphPrune::RemoveOrphantVertices() from AssemblyGraph
    void RemoveOrphanVertices(void)  {
      //GraphPrune gp;
      gp.RemoveOrphanVertices(graph_ptr_);
    }

    // wrapper function to call GraphPrune::ResolveOrientation() from AssemblyGraph
    void ResolveOrientation(void) {
      //GraphPrune gp;
      gp.ResolveOrientation(graph_ptr_);
    }

    // wrapper function to call GraphPrune::ResolveSequence() from AssemblyGraph
    void ResolveSequence(void)  {
      gp.ResoveSequence(graph_ptr_);
    }

  private:
    GraphPrune gp;
    //GraphTraversal gt;
};

#endif  //__ASSEMBLYGRAPH_H_
