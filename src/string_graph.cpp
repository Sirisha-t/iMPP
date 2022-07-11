#include "SG.h"
bool StrGraph::readFGSPredictions(char *fastaname, char *gffname){
  std::ifstream fin1, fin2;
  fin1.open(fastaname);
  clock_t t;
  if (!fin1.is_open())	    // fail to open
  {
    std::cout << "Error: " << fastaname << " doesnot exist!\n";
    return false;
  }
  else				// succeed to open
  {
    std::string line;
    int discard = 0;

    //boost::char_separator<char> sep(" \t");
    //boost::tokenizer<boost::char_separator<char> >::iterator it;
    /** timer **/
    t = clock();
    /***********/
    while(std::getline(fin1, line))
    {
      if(line[0] == '>')
      {
        line.erase(0,1);
        std::string read = line;
          if(readNameMap.count(read) == 0)
                readNameMap[read] = 0;
      }
      else
        discard++;
    }

  }
  t = clock() - t;
  std::cout << "Completed reading input fasta file in " << ((double)t)/CLOCKS_PER_SEC << " s\n";

  std::cout<<gffname<<"\n";
  fin2.open(gffname);

  if (!fin2.is_open())	    // fail to open
  {
    std::cout << "Error: " << gffname << " does not exist!\n";
    return false;
  }
  else				// succeed to open
  {
    std::string line;
    int discard = 0;
    int count=0;
    boost::char_separator<char> sep(" \t");
    boost::tokenizer<boost::char_separator<char> >::iterator it;
    /** timer **/

    t = clock();
    /***********/
    while(std::getline(fin2, line))
    {
      if(line[0] != '#')
      {
        boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
        it = tokens.begin();
        std::string read = *it;
        if(readNameMap[read] == 0)
        {
          readNameMap[read] = 1;
          count++;
        }
        //std::cout<<readNameMap[read]<<"\t"<<read<<"\n";

      }

    }
    t = clock() - t;
    std::cout << "Completed reading input gff file in " << ((double)t)/CLOCKS_PER_SEC << " s\n";
    std::cout<<"Pred reads count: "<<count<<"\n";
  }
  fin1.clear();
  fin1.close();
  fin2.clear();
  fin2.close();
  return true;
}

/***************************** for overlap graph ******************************/
bool StrGraph::readAsqgFile(char* filename)
{
  std::ifstream fin;
  StringVector p_header;
  String2Integer p_readID;
  std::string newname = filename;
      std::string outname = newname.substr(0,newname.find_last_of('.')) + ".wrong.asqg";
      std::ofstream fout(outname.c_str());

  // fin.open(config.IDIR_GRAPH.c_str());	// open file
  fin.open(filename);	// open file
  if (!fin.is_open())	    // fail to open
  {
    std::cout << "Error: " << filename << " doesnot exist!\n";
    return false;
  }
  else				// succeed to open
  {
    std::string line;
    boost::char_separator<char> sep(" \t");
    boost::tokenizer<boost::char_separator<char> >::iterator it;
    /** timer **/
    clock_t t;
    t = clock();
    /***********/
    while(std::getline(fin, line))
    {
      if(line[0] == 'V')
      {
        boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
        it = tokens.begin();            // cont[0], VT
        ++it;                           // cont[1], read name
        int id = p_readID.size();
        p_readID[*it] = id;
        predMap[id] = readNameMap[*it];
        //std::cout<<"asqg "<<*it<<"\t"<<p_asqgID[*it]<<"\n";
        ++it;                           // cont[2], sequence
        p_seq.push_back(*it);
      }
      if(line[0] == 'E')      break;
    }
    /** timer **/
    t = clock() - t;
    std::cout << "VT takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
    /***********/
    p_order = p_seq.size();
    /** remove double node storage **/
    p_traversed = BooleanVector(p_order, false);
    p_crossed   = BooleanVector(p_order, false);
    /** timer **/
    t = clock() - t;
    std::cout << "resize takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
    t = clock();
    /***********/
    do
    {
      if(line[0] == 'E')  // edge, sperate by '\t' and ' '
      {
        boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
        it = tokens.begin();            // cont[0], ED
        ++it;
        //std::string source = (*it);                         // cont[1], read1 name
        int source_id = p_readID[*it];
        ++it;                                                 // cont[2], read2 name
        int target_id = p_readID[*it];
        ++it;   
        if(std::stol(*it) > 100){
          fout<<line<<"\n";
        }                        // cont[3], A
        int a = std::stoi(*it);
        ++it;                           // cont[4], B
        int b = std::stoi(*it);
        ++it;                           // cont[5], l1
        int l1 = std::stoi(*it);
        ++it;                           // cont[6], C
        int c = std::stoi(*it);
        ++it;                           // cont[7], D
        int d = std::stoi(*it);
        ++it;                           // cont[8], l2
        int l2 = std::stoi(*it);
        ++it;                           // cont[9], "1" == RC
        bool ff_RC = *it == "1";

        /* add or locate the two vertices in the graph */
        BoostSTRVertex v_source, v_RC_source, v_target, v_RC_target;
        //std::cout<<source_id<<"\t"<<target_id<<"\n";

        if( predMap[source_id]!=1 || predMap[target_id]!=1)
        {
          //std::cout<<source<<"\t"<<readNameMap[source]<<"\t"<<source_id<<"\n";
        if(vertex.count(source_id) == 0)
        {
          STRVertexType node(source_id);
          node.len_ = l1;
          v_source = boost::add_vertex(node, *p_graph_);
          vertex[source_id] = v_source;
        }
        else
        {
          v_source = vertex[source_id];
        }

        /*if(vertex.count(source_id + p_order) == 0)
        {
          STRVertexType node(source_id + p_order);
          node.len_ = l1;
          v_RC_source = boost::add_vertex(node, *p_graph_);
          vertex[source_id + p_order] = v_RC_source;
        }
        else
        {
          v_RC_source = vertex[source_id + p_order];
        }*/

        if(vertex.count(target_id) == 0)
        {
          STRVertexType node(target_id);
          node.len_ = l2;
          v_target = boost::add_vertex(node, *p_graph_);
          vertex[target_id] = v_target;
        }
        else
        {
          v_target = vertex[target_id];
        }

        /*if(vertex.count(target_id + p_order) == 0)
        {
          STRVertexType node(target_id + p_order);
          node.len_ = l2;
          v_RC_target = boost::add_vertex(node, *p_graph_);
          vertex[target_id + p_order] = v_RC_target;
        }
        else
        {
          v_RC_target = vertex[target_id + p_order];
        }*/


       if(source_id != target_id)
        {
          /* add edge */
          if(ff_RC)         // r2 is reversed
          {
            if(a != 0)      // r1->r2
            {
              /* r1 -> r2' */
              std::pair<BoostSTREdge, bool> e_search = boost::add_edge(v_source, v_RC_target, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = a;
              }
              else

              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
              /* r2 -> r1' */
             e_search = boost::add_edge(v_target, v_RC_source, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = c;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
            }
            else            // r2->r1
            {
              /* r2' -> r1 */
              std::pair<BoostSTREdge, bool> e_search = boost::add_edge(v_RC_target, v_source, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = l2-1-d;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
              /* r1' -> r2 */
              e_search = boost::add_edge(v_RC_source, v_target, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = l1-1-b;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
            }
          }
          else
          {
            if(a != 0)     // r1->r2
            {
              /* r1 -> r2 */
             std::pair<BoostSTREdge, bool> e_search = boost::add_edge(v_source, v_target, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = a;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
              /* r2' -> r1' */
             /* e_search = boost::add_edge(v_RC_target, v_RC_source, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = l2-1-d;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }*/
            }
            else            // r2->r1
            {
              /* r2 -> r1 */
              std::pair<BoostSTREdge, bool> e_search = boost::add_edge(v_target, v_source, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = c;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
              /* r1' -> r2' */
             /* e_search = boost::add_edge(v_RC_source, v_RC_target, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = l1-1-b;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }*/
            }
          }
        }

    }
   
    }
  }
    while(std::getline(fin, line));
    /** timer **/
    t = clock() - t;
    std::cout << "ED takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
    /***********/

  }
  fin.clear();
  fin.close();
  return true;
}

void StrGraph::CondenseGraph()
{
  /** timer **/
  clock_t t;
  /***********/
  /* find seed vertex */
  t = clock();
  int total_node = 0;
  for(auto it = boost::vertices(*p_graph_).first; it != boost::vertices(*p_graph_).second; ++it)
  {
    //std::cout<<(*it)<<"\n";
    ++total_node;
    if(boost::in_degree(*it, *p_graph_) <= 0)
    { // orphant read is not seed
      if(boost::out_degree(*it, *p_graph_) > 0)   (*p_graph_)[*it].ff_seed = true;
    }
    else
    {
      if(boost::out_degree(*it, *p_graph_) != 1 )  (*p_graph_)[*it].ff_seed = true;
      if(boost::in_degree(*it, *p_graph_) > 1)    (*p_graph_)[*it].ff_seed = true;
    }
  }
  t = clock() - t;
  std::cout << total_node << " vertices in total\n";
  std::cout << "find seed takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";


  /* for each component */
  std::list<cycle_s> to_add_cycle;
  int cnt_cycle_node = 0;
  double condeseSE = 0, condenseCycle = 0, findComponent = 0;
  // remove double nodes -- test cz asqg
  BooleanVector connected = BooleanVector(p_order, false);
  for(auto it = boost::vertices(*p_graph_).first; it != boost::vertices(*p_graph_).second; ++it)
  {
    if(! p_crossed[ (*p_graph_)[*it].rid_ ]) // untraversed vertex, new component
    {
      /* find connected component and source_edges */
      t = clock();
      std::stack<BoostSTRVertex> to_visit;
      to_visit.push(*it);
      std::list<BoostSTRVertex> this_component;
      int cnt_vertex = 0;
      while(!to_visit.empty())
      {
        
        BoostSTRVertex top_vertex = to_visit.top();
        to_visit.pop();
        int rid = (*p_graph_)[top_vertex].rid_;
        //std::cout<<"top "<<rid<<"\t"<<connected[rid]<<"\n";

        if(!connected[rid])
        {
          connected[rid] = true;
          this_component.push_back(top_vertex);
          ++cnt_vertex;
          p_crossed[rid] = true;
          //if(rid < p_order) p_crossed[rid+p_order] = true;
          //else              p_crossed[rid-p_order] = true;
          for(auto it_v = boost::adjacent_vertices(top_vertex, *p_graph_).first; it_v != boost::adjacent_vertices(top_vertex, *p_graph_).second; ++it_v)
          {
            if(!connected[ (*p_graph_)[*it_v].rid_ ]) to_visit.push(*it_v);
          }

          for(auto it_v = boost::inv_adjacent_vertices(top_vertex, *p_graph_).first; it_v != boost::inv_adjacent_vertices(top_vertex, *p_graph_).second; ++it_v)
          {
            if(!connected[ (*p_graph_)[*it_v].rid_ ]) to_visit.push(*it_v);
          }
        }
      }
      t = clock() - t;
      findComponent += ((double)t)/CLOCKS_PER_SEC;

      /* condense this component */
      t = clock();
      bool ff_singleCycle = true;
      for(auto it_v = this_component.begin(); it_v != this_component.end(); ++it_v)
      {
        if((*p_graph_)[*it_v].IsSeed() && (!p_traversed[ (*p_graph_)[*it_v].rid_ ]))
        {
          ff_singleCycle = false;
          condense(*it_v, to_add_cycle);
        }
      }
      t = clock() - t;
      condeseSE += ((double)t)/CLOCKS_PER_SEC;

      t = clock();
      if(ff_singleCycle )  // if this component is single loop
      {
        cnt_cycle_node += cnt_vertex;
        auto it_v = this_component.begin();
        (*p_graph_)[*it_v].ff_seed = true;
        condense(*it_v, to_add_cycle);
      }
      t = clock() - t;
      condenseCycle += ((double)t)/CLOCKS_PER_SEC;
    }
  }
  std::cout << "Finding components took " << findComponent << " s\n";
  std::cout << "Consensing SE took " << condeseSE << " s\n";
  std::cout << "Condensing cycle took " << condenseCycle << " s\n";
  std::cout << cnt_cycle_node << " cycle nodes\n";
  /* modify graph */
  t = clock();
  std::list<BoostSTREdge> to_delete_edge;
  for(auto it_e = boost::edges(*p_graph_).first; it_e != boost::edges(*p_graph_).second; ++it_e)
  {
    if( (*p_graph_)[*it_e].ff_delete ) to_delete_edge.push_back(*it_e);
  }
  std::list<BoostSTRVertex> to_delete_vertex;
  for(auto it_v = boost::vertices(*p_graph_).first; it_v != boost::vertices(*p_graph_).second; ++it_v)
  {
    if( (*p_graph_)[*it_v].ff_delete ) to_delete_vertex.push_back(*it_v);
  }
  for(auto it = to_delete_edge.begin(); it != to_delete_edge.end(); ++it)     boost::remove_edge(*it, *p_graph_);
  for(auto it = to_delete_vertex.begin(); it != to_delete_vertex.end(); ++it) boost::remove_vertex(*it, *p_graph_);

  p_node_id = 0;
  for(auto it = to_add_cycle.begin(); it != to_add_cycle.end(); ++it)
  {
    BoostSTRVertex s = it->head;
    BoostSTRVertex t = it->tail;
    std::pair<BoostSTREdge, bool> edge_new = add_edge(s, t, *p_graph_);
    if(edge_new.second)
    {
      (*p_graph_)[edge_new.first].SetCondensedTag(true);
      (*p_graph_)[edge_new.first].path_info_ = it->path;
      (*p_graph_)[edge_new.first].ff_cycle   = it->ff_cycle;
      (*p_graph_)[edge_new.first].sid_       = p_node_id;
      ++p_node_id;
    }
  }

  t = clock() - t;
  std::cout << "Graph modification took " << ((double)t)/CLOCKS_PER_SEC << " s\n";
}

void StrGraph::condense(const BoostSTRVertex seed,
  std::list<cycle_s>& to_add_cycle)
  {
    for(auto it_e = boost::out_edges(seed, *p_graph_).first; it_e != boost::out_edges(seed, *p_graph_).second; ++it_e)
    condense(*it_e, to_add_cycle);
  }

  void StrGraph::condense(const BoostSTREdge source_edge,
    std::list<cycle_s>& to_add_cycle)
    {
      BoostSTREdge current_edge = source_edge;
      BoostSTRVertex head = boost::source(current_edge, *p_graph_);
      BoostSTRVertex tail = head;
      IntegerVector path_info;
      path_info.push_back((*p_graph_)[head].rid_);
      path_info.push_back((*p_graph_)[head].len_);
      while(1)
      {
        p_traversed[ (*p_graph_)[tail].rid_ ] = true;
        (*p_graph_)[current_edge].ff_delete = true;             //boost::remove_edge(current_edge, *p_graph_);
        if(tail != head)  (*p_graph_)[tail].ff_delete = true;   // boost::remove_vertex(to_delete, *p_graph_);

        tail = boost::target(current_edge, *p_graph_);          // define the new tail vertex
        path_info.push_back((*p_graph_)[current_edge].p_A);
        path_info.push_back((*p_graph_)[tail].rid_);
        path_info.push_back((*p_graph_)[tail].len_);

        if((*p_graph_)[tail].IsSeed())  break;
        else                            current_edge = *(boost::out_edges(tail, *p_graph_)).first;
      }
      path_info.push_back(-1);
      if( boost::out_degree(tail, *p_graph_) == 0)  p_traversed[ (*p_graph_)[tail].rid_ ] = true;

      cycle_s new_cycle;
      new_cycle.head = head;
      new_cycle.tail = tail;
      new_cycle.path = path_info;
      new_cycle.ff_cycle = ((*p_graph_)[head].rid_ == (*p_graph_)[tail].rid_ );
      to_add_cycle.push_back(new_cycle);
    }

    void StrGraph::writeGraph(char* filename)
    {
      /** timer **/
      clock_t t;
      t = clock();
      int cnt = 0;
      /***********/
      std::string newname = filename;
      std::string outname = newname.substr(0,newname.find_last_of('.')) + ".StringGraph.fq";
      std::ofstream fout(outname.c_str());
      for(auto it_e = boost::edges(*p_graph_).first; it_e != boost::edges(*p_graph_).second; ++it_e)
      {
        if((*p_graph_)[*it_e].IsCondensed())
        {
          ++cnt;
          std::string sequence = "";
          IntegerVector path_info = (*p_graph_)[*it_e].path_info_;
          fout << ">" << (*p_graph_)[*it_e].sid_ << ",";
          int p_before = 1;
          for(int i = 0; i < path_info.size(); i += 3)
          {
            int r_id = path_info[i];
            int r_len = path_info[i+1];
            int a = path_info[i+2];
            std::string temp = (r_id < p_order) ? p_seq[r_id] : RC(r_id-p_order);
            if(a > 0)   sequence += temp.substr(0, a);
            else        sequence += temp;

            fout << ((r_id < p_order) ? r_id : r_id-p_order) << ","
            << p_before << ","
            << r_len;
            if(i + 3 < path_info.size())  fout << ",";
            p_before = a + p_before;

            if(i + 6 >= path_info.size() && (*p_graph_)[*it_e].IsCycle())
            {
              i = -3;
              (*p_graph_)[*it_e].ff_cycle = false;
            }
          }

          BoostSTRVertex t = boost::target(*it_e, *p_graph_);
          auto it_out_e = boost::out_edges(t, *p_graph_).first;
          if(it_out_e != boost::out_edges(t, *p_graph_).second)
          {
            fout << ":" << (*p_graph_)[*it_out_e].sid_;
            ++it_out_e;
            while(it_out_e != boost::out_edges(t, *p_graph_).second)
            {
              fout << "," << (*p_graph_)[*it_out_e].sid_;
              ++it_out_e;
            }
          }
          fout << "\n";
          fout << sequence << "\n";
        }
      }
      /** timer **/
      t = clock() - t;
      std::cout << "write "<< cnt << " edges " << ((double)t)/CLOCKS_PER_SEC << " s\n";
      t = clock();
      cnt = 0;
      /***********/

      /* take care of orphant reads */
      for(int i = 0; i < p_order; ++i)
      {
        if(!p_traversed[i] && !p_traversed[i+p_order])
        {
          ++cnt;
          fout << ">" << p_node_id << "," << i << ",1," << p_seq[i].size() << "\n";
          fout << p_seq[i] << "\n";
          ++p_node_id;
        }
      }

      /** timer **/
      t = clock() - t;
      std::cout << "write "<< cnt << " orphants " << ((double)t)/CLOCKS_PER_SEC << " s\n";
      /***********/
      fout.close();
    }

    std::string StrGraph::RC(int i)
    {
      std::string res="";
      for(int j=0;j<p_seq[i].size();++j)
      {
        if(p_seq[i][j]=='A' || p_seq[i][j]=='a')
        res = "T" + res;
        else if (p_seq[i][j]=='T' || p_seq[i][j]=='t')
        res = "A" + res;
        else if (p_seq[i][j]=='G' || p_seq[i][j]=='g')
        res = "C" + res;
        else if (p_seq[i][j]=='C' || p_seq[i][j]=='c')
        res = "G" + res;
      }
      return res;
    }

    StrGraph::StrGraph()
    {
      p_graph_ = new BoostSTRGraph();
    }

    StrGraph::~StrGraph()
    {
      delete p_graph_;
    }
