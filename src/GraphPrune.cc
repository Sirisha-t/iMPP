//  Author: Cuncong Zhong
//  Last modification: 11/10/2021

#include "GraphPrune.h"

using namespace std;

//==========================================================================================
//      Public methods implementations
//==========================================================================================

// remove orphant vertices
// parameter:
//    g: the pointer to the graph where the function operate on
void GraphPrune::RemoveOrphanVertices(AssemblyGraphType *g)    {
    // remove orphan vertices
    auto it_v = boost::vertices(*g).first;
    while(it_v != boost::vertices(*g).second) {
        BoostNodeType v = *it_v; ++ it_v;
        if(boost::in_degree(v, *g) <= 0 && boost::out_degree(v, *g) <= 0) {
            boost::remove_vertex(v, *g);
        }
    }
    return;
}

// Condence single paths
// parameter:
//    g: the pointer to the graph where the function operates on
void GraphPrune::CondenceSinglePaths(AssemblyGraphType *g)  {

    // TODO: the following code has not been tested
    /*
    std::list<BoostSTREdge> source_edges;
    auto it = boost::vertices(*p_graph_).first;
    while(it != boost::vertices(*p_graph_).second) {
        if(boost::in_degree(*it, *p_graph_) <= 0) {
            auto it_e = boost::out_edges(*it, *p_graph_).first;
            while(it_e != boost::out_edges(*it, *p_graph_).second) {
                source_edges.push_back(*it_e); ++ it_e;
            }
        }
        ++ it;
    }
  
    // condense the graph
    for(auto it = source_edges.begin(); it != source_edges.end(); ++ it) {
        Condense(seq, *it);
    }
  
    // double check if all edges are filled out
    // if not, the edge is a cycle, break it
    auto it_e = boost::edges(*p_graph_).first;
    while(it_e != boost::edges(*p_graph_).second) {
        bool to_del = false; BoostSTREdge de;
        if((*p_graph_)[*it_e].path_info_.size() <= 0)  {
            de = *it_e; to_del = true;
        }
        ++ it_e;
        if(to_del)  boost::remove_edge(de, *p_graph_);
    }
    return;
    */
}

// Resolve read orientation; for each read only one direction is retained in the assembly graph
// Algorithm: Start with an unresolved node with the highest degree and set it to positive strand
// Repeat:  set its neigbhors' orientation according to the edge;
//          if the neighbor has already been set and is conflicting with the current edge, remove the current edge
//          mark the neighbors as resolved
//          return to step 1
// Internal data structure used: priority_queue
void GraphPrune::ResolveOrientation(AssemblyGraphType *g)   {

    // TODO: address ov_olp_ values when resolving orientations

    // defining the comparison class for priority queue
    // each element is a pair <BoostNodeType, int>;
    // the first field is the node, the second field is the degree of the node
    class cmp_node_degree   {
      public:  
        // comparison operator
        bool operator() (const std::pair<BoostNodeType, int> &lhs, const std::pair<BoostNodeType,int> &rhs) const {
            // implement a MAX queue; assign higher priority to nodes with higher degree
            return lhs.second < rhs.second ? true : false;
        } 
    };
    // insert all nodes into priority queue
    std::priority_queue<std::pair<BoostNodeType, int>, std::vector<std::pair<BoostNodeType, int> >,  cmp_node_degree> degree_rank;
    auto it_v = boost::vertices(*g).first;
    while(it_v != boost::vertices(*g).second) {
        degree_rank.push(std::make_pair(*it_v, boost::out_degree(*it_v, *g)));
        ++ it_v;
    }

    // DEBUG
    //IDType g_ck_source = 286042;
    //IDType g_ck_target = 296777;

    // process the nodes one-by-one 

    std::vector<BoostEdgeType> to_delete;                           // the set of edges to be deleted
    std::vector<std::pair<BoostNodeType, BoostNodeType> > to_add;   // the set of edges to be added
    std::vector<GraphEdgeType> to_add_info;                         // the content of the edge to be added  
    size_t count =0;
    while(!degree_rank.empty()) {

        count++;
        std::pair<BoostNodeType, int> node_info = degree_rank.top();
        degree_rank.pop();
        //cout<<"DEBUG: node id: "<<(*g)[node_info.first].id_<<endl;
        //cout<<"DEBUG: Loop count: "<<count<<endl;
        //cout << "DEBUG: node degree:    " << boost::out_degree(node_info.first, *g) << endl;
        if((*g)[node_info.first].IsResolved()) continue;    // skip the source if it is resolved
        // set positive orientation for the source
        (*g)[node_info.first].SetOrientation(true);
        (*g)[node_info.first].SetResolved(true);
        // use BFS style traversal to resolve orientation
        std::queue<BoostNodeType> n_set;    // node set
        n_set.push(node_info.first);
        while(!n_set.empty()) {
            BoostNodeType s = n_set.front();  n_set.pop();
            // making sure each node can be used as a source for only once
            // otherwise may cause infinite loop if cycles exist
            if((*g)[s].IsVisited())    {
                continue;
            }   else    {
                (*g)[s].SetVisited(true); 
            }
            bool s_ori = (*g)[s].GetOrientation();   
            

              // DEBUG
            //cout << "DEBUG: ******" << endl;
            //cout << "DEBUG: source resolved?:   " << (*g)[s].IsResolved() << endl;
            //cout << "DEBUG: source orientation:   " << (*g)[s].GetOrientation() << endl;
            //(*g)[s].PrintInfo();
            //cout << "DEBUG: ******" << endl;
          /*  if((*g)[s].id_ == g_ck_source)    {
                cout << "DEBUG: " << g_ck_source << " used as a source node" << endl;
                cout << "DEBUG: " << g_ck_source << " orientation:    " << (*g)[s].orientation_ << "  is visited: " << (*g)[s].visited_ << "   is resolved:    " << (*g)[s].resolved_ << endl;
            }   else if ((*g)[s].id_ == g_ck_target)  {
                cout << "DEBUG: " << g_ck_target << " used as a source node" << endl;
                cout << "DEBUG: " << g_ck_target << " orientation:    " << (*g)[s].orientation_ << "  is visited: " << (*g)[s].visited_ << "   is resolved:    " << (*g)[s].resolved_ << endl;
            }
            */

              
            // iterate through all edges
            auto it_e = boost::out_edges(s, *g).first;
            while(it_e != boost::out_edges(s, *g).second) {

                BoostEdgeType cr_e = *it_e; ++ it_e;

                if((*g)[cr_e].IsVisited())    {
                    continue;       // if the edge has been visited, skip the edge
                }
                (*g)[cr_e].SetVisited(true);        // marking the edge as visited
                BoostNodeType t = boost::target(cr_e, *g);

                  // DEBUG
                //cout << "DEBUG: ******" << endl;
                //cout<<"DEBUG: target id: "<<(*g)[t].id_<<endl();
                //cout << "DEBUG: target resolved?:   " << (*g)[t].IsResolved() << endl;
                //cout << "DEBUG: target orientation:   " << (*g)[t].GetOrientation() << endl;
                //(*g)[t].PrintInfo();
                //cout << "DEBUG: ******" << endl;
                /*if((*g)[t].id_ == g_ck_source)    {
                    cout << "DEBUG: " << g_ck_source << " used as a target node" << endl;
                    cout << "DEBUG: " << g_ck_source << " orientation:    " << (*g)[t].orientation_ << "  is visited: " << (*g)[t].visited_ << "   is resolved:    " << (*g)[t].resolved_ << endl;
                }   else if ((*g)[t].id_ == g_ck_target)  {
                    cout << "DEBUG: " << g_ck_target << " used as a target node" << endl;
                    cout << "DEBUG: " << g_ck_target << " orientation:    " << (*g)[t].orientation_ << "  is visited: " << (*g)[t].visited_ << "   is resolved:    " << (*g)[t].resolved_ << endl;
                }*/
                

                // if target not resolved, assign orientation   
                bool is_e_removed = false;      // a tag indicating whether the edge has been removed
                if(!(*g)[t].IsResolved())   {
                    //cout << "DEBUG: target not resolved, assigning orientation." << endl;
                    (*g)[t].SetOrientation((*g)[cr_e].IsRevComplement() ? !s_ori : s_ori);     // setting orientation
                    (*g)[t].SetResolved(true);      // marking the target node as resolved
                    n_set.push(t);                  // adding the target node to the queue
                }
                // if target resolved  
                else {
                    bool t_ori = (*g)[t].GetOrientation();
                    if((s_ori == t_ori) == !(*g)[cr_e].IsRevComplement())  {   // nodes orientations are consistent with edge info
                        //cout << "DEBUG: target resolved and consistent, do nothing." << endl;
                        ;    // if consistent, do nothing
                    }   else    {
                        // if not consistent, mark edge for deletion
                        //cout << "DEBUG: target resolved and inconsistent, deleting edge." << endl;
                        to_delete.push_back(cr_e);
                        is_e_removed = true;
                    }
                }

                // handle the overlapped region: flip if the orientation is negative strand
                if(!is_e_removed)   {
                    if(!s_ori)    {     // negative strand
                        SeqIdxType len = (*g)[s].GetSeqLen();
                        SeqIdxType b = (*g)[cr_e].ov_pos_[0];
                        SeqIdxType e = (*g)[cr_e].ov_pos_[1];
                        (*g)[cr_e].ov_pos_[0] = len - e - 1;
                        (*g)[cr_e].ov_pos_[1] = len - b - 1;
                    }
                    if(!(*g)[t].GetOrientation())    {      // negative strand
                        SeqIdxType len = (*g)[t].GetSeqLen();
                        SeqIdxType b = (*g)[cr_e].ov_pos_[2];
                        SeqIdxType e = (*g)[cr_e].ov_pos_[3];
                        (*g)[cr_e].ov_pos_[2] = len - e - 1;
                        (*g)[cr_e].ov_pos_[3] = len - b - 1;
                    }
                }

                // if the source strand is different from the edge information and the edge is not removed, needs to flip the edge direction
                if(s_ori != (*g)[cr_e].GetSrcOrientation() && !is_e_removed)  {
                    to_delete.push_back(cr_e);                  // delete the current edge
                    to_add.push_back(std::make_pair(t, s));     // add the reverse direction
                    // swap the overlap position information
                    SeqIdxType tmp1 = (*g)[cr_e].ov_pos_[0];
                    SeqIdxType tmp2 = (*g)[cr_e].ov_pos_[1];
                    (*g)[cr_e].ov_pos_[0] = (*g)[cr_e].ov_pos_[2];
                    (*g)[cr_e].ov_pos_[1] = (*g)[cr_e].ov_pos_[3];
                    (*g)[cr_e].ov_pos_[2] = tmp1;
                    (*g)[cr_e].ov_pos_[3] = tmp2;
                    to_add_info.push_back((*g)[cr_e]);          // record the information

                    // TODO: need to flip CIGAR information

                      // DEBUG
                 /*   GraphEdgeType tmp_e = (*g)[cr_e];       // record the edge information
                    boost::remove_edge(cr_e, *g);           // remove the existing edge
                    pair<BoostEdgeType, bool> e_add = boost::add_edge(t, s, *g);    // add an edge with reverse direction
                    if(!e_add.second)    {
                        cout << "Error: GraphPrune::ResolveOrientation: Failed to reverse edge direction! Abort." << endl;
                        exit(1);
                    }
                    (*g)[e_add.first] =  tmp_e; 
                    if((*g)[s].id_ == g_ck_source && (*g)[t].id_ == g_ck_target)    {
                        cout << "DEBUG: recording wrong edge:   " << (*g)[t].id_ << "  " << (*g)[s].id_ << endl; 
                    }
                    if((*g)[s].id_ == g_ck_target && (*g)[t].id_ == g_ck_source)    {
                        cout << "DEBUG: recording correct edge: " << (*g)[t].id_ << "  " << (*g)[s].id_ << endl; 
                    }
                    */
                    
                    
                }
            }
        }
    }
    cout<<"DEBUG: degree_rank loop count : "<<count<<endl;
    cout << "DEBUG: delete size:    " << to_delete.size() << endl;
    cout << "DEBUG: add size:   " << to_add.size() << endl;

    // delete all the edges
    for(size_t i = 0; i < to_delete.size(); ++ i)   {
        
        /*  // DEBUG
        BoostNodeType ck_s = boost::source(to_delete[i], *g);
        BoostNodeType ck_t = boost::target(to_delete[i], *g);
        if((*g)[ck_s].id_ == g_ck_source && (*g)[ck_t].id_ == g_ck_target)    {
            cout << "DEBUG: removing correct edge." << endl; 
        }   else if((*g)[ck_s].id_ == g_ck_target && (*g)[ck_t].id_ == g_ck_source)    {
            cout << "DEBUG: removing wrong edge." << endl; 
        }
        */

        boost::remove_edge(to_delete[i], *g);                                   // delete the edge from the graph

    }
    // add all the edges
    assert(to_add.size() == to_add_info.size());
    for(size_t i = 0; i < to_add.size(); ++ i)   {

        /*  // DEBUG
        BoostNodeType ck_s = to_add[i].first;
        BoostNodeType ck_t = to_add[i].second;
        if((*g)[ck_s].id_ == g_ck_source && (*g)[ck_t].id_ == g_ck_target)    {
            cout << "DEBUG: adding correct edge." << endl; 
        }   else if((*g)[ck_s].id_ == g_ck_target && (*g)[ck_t].id_ == g_ck_source)    {
            cout << "DEBUG: adding wrong edge." << endl; 
        }
        */

        to_add_info[i].SetSrcOrientation((*g)[to_add[i].first].GetOrientation());   // update the source orientation information
        boost::add_edge(to_add[i].first, to_add[i].second, to_add_info[i], *g);     // add the edge into the graph
    }
    return;
}

// Resolve the sequences in the assembly graph.
// Use reverse complementary if the read is assigned negative strand
// Also update all edges as no reverse-completementary
// Parameter:
//    g: the pointer to the graph where the function operates on
void GraphPrune::ResoveSequence(AssemblyGraphType *g)   {
    // check the consisency of the graph
    auto it_e = boost::edges(*g).first;
    while(it_e != boost::edges(*g).second) {

        /*
        BoostEdgeType ck_e = *it_e;
        BoostNodeType ck_s = boost::source(ck_e, *g);
        BoostNodeType ck_t = boost::target(ck_e, *g);
        if((*g)[ck_s].id_ == 296777 && (*g)[ck_t].id_ == 286042)    {
            cout << "DEBUG: wrong edge direction before ResolveSequence!!!" << endl;
        }   else if ((*g)[ck_t].id_ == 296777 && (*g)[ck_s].id_ == 286042)  {
            cout << "DEBUG: correct edge direction before ResolveSequence!!!" << endl;
        }
        */


       if(((*g)[*it_e].is_rc_ && (*g)[boost::source(*it_e, *g)].orientation_ != (*g)[boost::target(*it_e, *g)].orientation_) ||
         !((*g)[*it_e].is_rc_ && (*g)[boost::source(*it_e, *g)].orientation_ == (*g)[boost::target(*it_e, *g)].orientation_)) {
            // consistent edge and nodes
            (*g)[*it_e].is_rc_ = false;
       }    else    {
           // inconsistency found
            cerr << "Error:  GraphPrune::ResolveSequence: Inconsistent edge/vertex orientation found. Call ResolveOrientation() before calling ResolveSequence(). Exit program." << endl;
            exit(1);        
       }
        ++ it_e;
    }
    // convert the sequence for each vertex 
    auto it_v = boost::vertices(*g).first;
    while(it_v != boost::vertices(*g).second) {
        if(!(*g)[*it_v].orientation_)    {  // if reverse complementary
            (*g)[*it_v].orientation_ = true;
            StringUtils::InplaceRevComp((*g)[*it_v].str_, (*g)[*it_v].str_len_);
        }
        ++ it_v;
    }
    return;
}
