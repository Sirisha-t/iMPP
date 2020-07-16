#include "SG.h"

bool StrGraph::loadAnchor(
  char* pEdgeSet,
  int MAX_EXTEND_LENGTH,
  anchors& anchor_set)
{
  std::ifstream fin;
  fin.open(pEdgeSet);
  if (!fin.is_open())	// fail to open
  {
    std::cout << "Error: unable to open " << pEdgeSet << "\n";
    return false;
  }
  else				// succeed to open
  {
    std::string line;
    IntegerVector pred_nodes;
    IntegerVector anc;
    //std::vector<int>::iterator an;
    boost::char_separator<char> sep("\t");
    boost::tokenizer<boost::char_separator<char> >::iterator it;
    while(std::getline(fin, line))
    {
      if(line[0] != '#')
      {
        //boost::split(fields,line, boost::is_any_of("\t"));
        //std::string edge = fields[0];
        boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
        it = tokens.begin();
        int node = std::stoi(*it);
        pred_nodes.push_back(node);
      }
    }
    IntegerSet p_nodes_set(p_all_nodes.begin(), p_all_nodes.end());
    IntegerSet pred_nodes_set(pred_nodes.begin(), pred_nodes.end());
    //std::sort(pred_nodes.begin(), pred_nodes.end());
    //std::sort(p_nodes_set.begin(), p_nodes_set.end());
    std::set_difference(p_nodes_set.begin(),p_nodes_set.end(),\
                        pred_nodes.begin(),pred_nodes.end(),\
                        std::back_inserter(anc));
    for(auto an = anc.begin(); an!= anc.end(); ++an)
    {
      struct anchor temp;
      parseAnchorLine((*an),temp, MAX_EXTEND_LENGTH);
      anchor_set.push_back(temp);
      //std::cout<<anc.size()<<"\t"<<*an<<"\n";
    }

    }

  fin.clear();
  fin.close();

  return true;
}

bool StrGraph::parseAnchorLine(
  int & node,
  struct anchor& new_anchor,
  int MAX_EXTEND_LENGTH)
{
    // query name   : 595,11855,1,2003,6220,1948,246:593,594
    //boost::char_separator<char> sep(","); //lines from here to int node- latest addition
  //  boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
  //  boost::tokenizer<boost::char_separator<char> >::iterator it;
    int length = p_seq[node].size();
    //std::cout<<"Length of node is: "<<length<<"\n";
    int left, right;
    //new_anchor.push_back(node);
    if(MAX_EXTEND_LENGTH > 0)
    {
        left = MAX_EXTEND_LENGTH;
        right = MAX_EXTEND_LENGTH;
    }
    new_anchor.path = IntegerList{node};
    new_anchor.M_upstream = left;
    new_anchor.M_downstream = right;

    //std::cout<<"Path in anchor is :"<<new_anchor.path<<"\n";

  return true;
}

bool StrGraph::loadPredictedEdge(
  char* gffname,
  predEdges& pEdge_set)
{

  std::ifstream fin;
  fin.open(gffname);

  if (!fin.is_open())	// fail to open
  {
    std::cout << "Error: unable to open " << gffname << "\n";
    return false;
  }
  else				// succeed to open
  {
      std::string line;
      while(std::getline(fin, line))
      {
        if(line[0] != '#')
        {
          struct pEdgeInfo temp;
          int pnode = parseGFFLine(line, temp);
          if(pEdge_set.count(pnode) == 0)
            pEdge_set[pnode] = temp;
          //anchor_set = temp;
        }
      }
  }
  fin.clear();
  fin.close();

  return true;
}

int StrGraph::parseGFFLine(
  const std::string& line,
  struct pEdgeInfo& new_edge)
{
    // query name   : 595,11855,1,2003,6220,1948,246:593,594
    boost::char_separator<char> sep("\t"); //lines from here to int node- latest addition
    boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
    boost::tokenizer<boost::char_separator<char> >::iterator it;
    it = tokens.begin();
    int node = std::stoi(*it);
    if(pEdge.count(node) == 0)
    {
      pEdge[node] = pEdge.size() + 1;
      new_edge.node_id = node;
      ++it;
      ++it;
      ++it;
      new_edge.start = std::stoi(*it);
      ++it;
      new_edge.end = std::stoi(*it);
      ++it;
      ++it;
      if((*it) == "+")
        new_edge.dir = 1; //1 for + strand
      else
        new_edge.dir = 0; // 0 for - strand
      ++it;
      new_edge.frame = std::stoi(*it);
    }
    //std::cout<<new_edge.node_id<<"\t"<<new_edge.start<<"\t"<<new_edge.end<<"\t"<<new_edge.dir<<"\t"<<new_edge.frame<<"\n";

  return node;
}

bool StrGraph::DirectedDFS(
  char* fgsoutname,
  int MAX_EXTEND_LENGTH,
  char* path_name)
{
  /* load Anchor */
  /** timer **/
  clock_t t;
  t = clock();
  /***********/
  anchors anchor_set;
  predEdges prededge_set;
  if(!loadAnchor(fgsoutname,MAX_EXTEND_LENGTH, anchor_set)) return false;
  /** timer **/
  t = clock() - t;
  std::cout << "Loading anchor takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  t = clock();
  /***********/
  if(!loadPredictedEdge(fgsoutname,prededge_set)) return false;
  /** timer **/
  t = clock() - t;
  std::cout << "Loading predicted edges takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  t = clock();

  /* DFS */
  // clear file
  std::ofstream path_out;
  path_out.open(path_name);
  path_out.clear();
  path_out.close();
  p_node_id = 0;
  p_crossed = BooleanVector(p_seq.size(), false);
  int anc_count=0;
  size_t total_path = 0;

  for(auto it = anchor_set.begin(); it != anchor_set.end(); ++it)
  {
    int anc_path=0;
    anc_count++;
    //std::cout<<"Anchor to be traversed: "<<anc_count<<"\n";
    DirectedDFS((*it),prededge_set, MAX_EXTEND_LENGTH, path_name, anc_path);
    size_t total_path = total_path + anc_path;
  }

  /** timer **/
  t = clock() - t;
  std::cout << "Extending anchor takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  /***********/
  return true;
}

void StrGraph::DirectedDFS(struct anchor& anchor, predEdges& prededge_set, int MAX_EXTEND_LENGTH, char* path_name, int& pathnum)
{
  std::stack<struct anchor> checking_stacks, half_done_stacks;
  checking_stacks.push(anchor);
  // go up first
  while(!checking_stacks.empty())
  {
    struct anchor top_anchor = checking_stacks.top();
    checking_stacks.pop();
    if(top_anchor.M_upstream > 0 ) // need to go deeper
    {
      IntegerList path = top_anchor.path;
      int first_node_id = path.front();
      BoostSTRVertex first_node = vertex[first_node_id];
      auto it_v_inv = boost::inv_adjacent_vertices(first_node, *p_graph_).first;
      if(it_v_inv == boost::inv_adjacent_vertices(first_node, *p_graph_).second)  // no head
      {
        half_done_stacks.push(top_anchor);
      }
      else                                                                        // for each head
      {
        for(; it_v_inv != boost::inv_adjacent_vertices(first_node, *p_graph_).second; ++it_v_inv)
        {
          int new_head = (*p_graph_)[*it_v_inv].rid_;
          //std::cout<<"Head: "<<new_head<<"\t";
          if(!formCycle(path, new_head))  // no cycle allowed
          {
            path.push_front(new_head);
            p_crossed[new_head] = true;
            int last_read_length;
            if(pEdge.count(new_head) == 0)
            {
                IntegerVector reads = p_read_on_node[new_head];
                int last_read_id = *(reads.end()-1);
                last_read_length = p_read_length[last_read_id];
              // update anchor stack
              struct anchor new_anchor;
              new_anchor.path = path;
              new_anchor.M_upstream = top_anchor.M_upstream - p_seq[new_head].size() + last_read_length;
              new_anchor.M_downstream = top_anchor.M_downstream;

              checking_stacks.push(new_anchor);
              path.pop_front();
            }
            else
            {
                  if((prededge_set[new_head].end -  prededge_set[new_head].start + 1 ) <= MAX_EXTEND_LENGTH)
                   {   /*struct anchor new_anchor;
                      new_anchor.path = path;
                      new_anchor.M_upstream = 0;
                      new_anchor.M_downstream = top_anchor.M_downstream;
                      checking_stacks.push(new_anchor);*/
                      top_anchor.path = path;
                      half_done_stacks.push(top_anchor);
                      path.pop_front();
                    }
            }
           //std::cout<<"Upstream value: "<<new_anchor.M_upstream<<"\n";
        }
      }
    }
  }
  else                          // no need to go deeper
  {
    half_done_stacks.push(top_anchor);
  }
}
  // go down
  while(!half_done_stacks.empty())
  {
    struct anchor top_anchor = half_done_stacks.top();
    half_done_stacks.pop();
      if(top_anchor.M_downstream > 0)
      {
        IntegerList path = top_anchor.path;
        int last_node_id = path.back();
        BoostSTRVertex last_node = vertex[last_node_id];
        auto it_v = boost::adjacent_vertices(last_node, *p_graph_).first;
        if(it_v == boost::adjacent_vertices(last_node, *p_graph_).second)  // no head
        {
                  pathnum++;
                  //std::cout<<"SAVE PATH\n";
                 savePath(path_name, top_anchor);
        }
        else
        {
          for(; it_v != boost::adjacent_vertices(last_node, *p_graph_).second; ++it_v)
          {
            int new_tail = (*p_graph_)[*it_v].rid_;
            //std::cout<<"Tail: "<<new_tail<<"\t";
            if(!formCycle(path, new_tail))  // no cycle allowed
            {
              path.push_back(new_tail);
              p_crossed[new_tail] = true;
              int first_read_length;
              if(pEdge.count(new_tail) == 0)
              {
                  IntegerVector reads = p_read_on_node[new_tail];
                  int first_read_id = reads[0];
                  first_read_length = p_read_length[first_read_id];
                // update anchor stack
                struct anchor new_anchor;
                new_anchor.path = path;
                new_anchor.M_upstream = top_anchor.M_upstream;
                new_anchor.M_downstream = top_anchor.M_downstream - p_seq[new_tail].size() + first_read_length;

                half_done_stacks.push(new_anchor);
                path.pop_back();
              }
              else{
                      if(prededge_set[new_tail].start <= MAX_EXTEND_LENGTH &&  prededge_set[new_tail].end <= MAX_EXTEND_LENGTH)
                        {
                          /*struct anchor new_anchor;
                          new_anchor.path = path;
                          new_anchor.M_upstream = top_anchor.M_upstream;
                          new_anchor.M_downstream = 0;*/

                          //half_done_stacks.push(new_anchor);
                          top_anchor.path = path;
                          savePath(path_name, top_anchor);
                          path.pop_back();
                        }
              }

            }
        }
      }
    }
    else{
          pathnum++;
          //std::cout<<"SAVE PATH \n";
          savePath(path_name, top_anchor);
    }
  }
 // std::cout<<"Number of paths per anchor: "<<pathnum<<"\n";
}


void StrGraph::savePath(char* path_name, const struct anchor& anchor)
{
  ++p_node_id;
  std::ofstream path_out;
  path_out.open(path_name,std::ios::in|std::ios::out|std::ios::app);
  IntegerList path = anchor.path;
  std::string seq;
  path_out << ">" << p_node_id;
  for(auto it = path.begin(); it != path.end(); ++it) path_out << "," << (*it);
  path_out << " ";
  // eg: path_i,node_x1,node_x2,...,node_xn r_y1,begin,length,r_y2,begin,length,..,r_ym,begin,length,
  // all integer, don't cut here
  int M = 0;
  for(auto it = path.begin(); it != path.end(); ++it)
  {
    int edge = *it;
    // assembly sequence
    std::string full = p_seq[edge];
    int first_read = p_read_on_node[edge][0];
    if(it != path.begin()) seq += full.substr(p_read_length[first_read]);
    else                   seq += full;
    // print reads on path
    IntegerVector reads = p_read_on_node[edge];
    IntegerVector locus = p_pos_on_node[edge];
    if(it != path.begin())
    {
      for(int i = 1; i < reads.size(); i++)
      {
        path_out  << reads[i]                     << ","
                  << std::to_string(locus[i] + M) << ","
                  << p_read_length[reads[i]]      << ",";
      }
    }
    else
    {
      for(int i = 0; i < reads.size(); i++)
      {
        path_out  << reads[i]                     << ","
                  << std::to_string(locus[i] + M) << ","
                  << p_read_length[reads[i]]      << ",";
      }
    }
    M += full.size() - p_read_length[*(reads.end()-1)];
  }
  path_out << "\n" << seq << "\n";
  path_out.close();

}


bool formCycle(IntegerList & path, int new_tail)
{
  for(auto it = path.begin(); it != path.end(); ++it)
  {
    if((*it) == new_tail)  return true;
  }
  return false;
}


bool StrGraph::readFqFile(char* filename)
{
  /** timer **/
  clock_t t;
  t = clock();
  /***********/
  std::ifstream fin;
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
    boost::char_separator<char> sep(",");
    boost::tokenizer<boost::char_separator<char> >::iterator it;
    while(std::getline(fin, line))
    {
      if(line[0] == '>')
      {
        line.erase(0,1);
        int colon = line.find(':');
        if(colon == std::string::npos)  // no tail
        {
          std::string header = line.substr(0,colon);
          int colon = line.find(',');
          std::string path_info_str = header.substr(colon+1);
          // p_header.push_back(path_info_str);

          boost::tokenizer<boost::char_separator<char> > tokens(header, sep);
          it = tokens.begin();
          int node_id = std::stoi(*it);
          p_all_nodes.push_back(node_id);
          ++it;
          int read_id = -1;
          IntegerVector read_on_this_node, loc_on_this_node;
          for(int i = 0 ; it != tokens.end(); ++it, ++i)
          /* format: read_id, location, read_length */
          {
            if(i % 3 == 0)        read_on_this_node.push_back(std::stoi(*it));
            else if( i % 3 == 1)  loc_on_this_node.push_back(std::stoi(*it));
            else                  p_read_length[*(read_on_this_node.end()-1)] = std::stoi(*it);
          }
          p_read_on_node.push_back(read_on_this_node);
          p_pos_on_node.push_back(loc_on_this_node);

          BoostSTRVertex v_source;
          if(vertex.count(node_id) == 0)
          {
            STRVertexType node(node_id);
            v_source = boost::add_vertex(node, *p_graph_);
            vertex[node_id] = v_source;
          }
          else
          {
            v_source = vertex[node_id];
          }
        }
        else
        {
          std::string header = line.substr(0,colon);
          std::string tail = line.substr(colon+1);
          int colon = line.find(',');
          std::string path_info_str = header.substr(colon+1);
          // p_header.push_back(path_info_str);

          boost::tokenizer<boost::char_separator<char> > tokens(header, sep);
          it = tokens.begin();
          int node_id = std::stoi(*it);
          p_all_nodes.push_back(node_id);
          ++it;
          IntegerVector read_on_this_node, loc_on_this_node;
          for(int i = 0 ; it != tokens.end(); ++it, ++i)
          /* format: read_id, location, read_length */
          {
            if(i % 3 == 0)        read_on_this_node.push_back(std::stoi(*it));
            else if( i % 3 == 1)  loc_on_this_node.push_back(std::stoi(*it));
            else                  p_read_length[*(read_on_this_node.end()-1)] = std::stoi(*it);
          }
          p_read_on_node.push_back(read_on_this_node);
          p_pos_on_node.push_back(loc_on_this_node);

          BoostSTRVertex v_source, v_target;
          if(vertex.count(node_id) == 0)
          {
            STRVertexType node(node_id);
            v_source = boost::add_vertex(node, *p_graph_);
            vertex[node_id] = v_source;
          }
          else
          {
            v_source = vertex[node_id];
          }

          boost::tokenizer<boost::char_separator<char> > tokens1(tail, sep);
          for(it = tokens1.begin(); it != tokens1.end(); ++it)
          {
            int target_id = std::stoi(*it);
            if(vertex.count(target_id) > 0)
            {
              v_target = vertex[target_id];
            }
            else
            {
              STRVertexType node(target_id);
              v_target = boost::add_vertex(node, *p_graph_);
              vertex[target_id] = v_target;
            }

            std::pair<BoostSTREdge, bool> e_search = boost::add_edge(v_source, v_target, *p_graph_);
            if(!e_search.second)          std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
          }
        }
      }
      else
      {
        p_seq.push_back(line);
      }
    }
  }
  p_order = p_seq.size();
  ff_stringGraph = true;
  fin.clear();
  fin.close();
  /** timer **/
  t = clock() - t;
  std::cout << "load fastq takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  /***********/
  return true;
}

StrGraph::StrGraph()
{
  p_graph_ = new BoostSTRGraph();
}

StrGraph::~StrGraph()
{
  delete p_graph_;
}

/*

bool StrGraph::readAsqgFile(char* filename)
{
  std::ifstream fin;
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
    clock_t t;
    t = clock();
    while(std::getline(fin, line))
    {
      if(line[0] == 'V')
      {
        boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
        it = tokens.begin();            // cont[0], VT
        ++it;                           // cont[1], read name
        int id = std::stoi(*it);
        ++it;                           // cont[2], sequence
        p_seq.push_back(*it);
      }
      if(line[0] == 'E')      break;
    }
    t = clock() - t;
    std::cout << "VT takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
    p_order = p_seq.size();
    p_traversed = BooleanVector(2*p_order, false);
    p_crossed   = BooleanVector(2*p_order, false);
    t = clock() - t;
    std::cout << "resize takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
    t = clock();
    do
    {
      if(line[0] == 'E')  // edge, sperate by '\t' and ' '
      {
        boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
        it = tokens.begin();            // cont[0], ED
        ++it;                           // cont[1], read1 name
        int source_id = std::stoi(*it);
        ++it;                           // cont[2], read2 name
        int target_id = std::stoi(*it);
        ++it;                           // cont[3], A
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

        BoostSTRVertex v_source, v_RC_source, v_target, v_RC_target;
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

        if(vertex.count(source_id + p_order) == 0)
        {
          STRVertexType node(source_id + p_order);
          node.len_ = l1;
          v_RC_source = boost::add_vertex(node, *p_graph_);
          vertex[source_id + p_order] = v_RC_source;
        }
        else
        {
          v_RC_source = vertex[source_id + p_order];
        }

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

        if(vertex.count(target_id + p_order) == 0)
        {
          STRVertexType node(target_id + p_order);
          node.len_ = l2;
          v_RC_target = boost::add_vertex(node, *p_graph_);
          vertex[target_id + p_order] = v_RC_target;
        }
        else
        {
          v_RC_target = vertex[target_id + p_order];
        }

        if(source_id != target_id)
        {
          if(ff_RC)         // r2 is reversed
          {
            if(a != 0)      // r1->r2
            {
              std::pair<BoostSTREdge, bool> e_search = boost::add_edge(v_source, v_RC_target, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = a;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
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
              std::pair<BoostSTREdge, bool> e_search = boost::add_edge(v_RC_target, v_source, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = l2-1-d;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
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
              std::pair<BoostSTREdge, bool> e_search = boost::add_edge(v_source, v_target, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = a;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
              e_search = boost::add_edge(v_RC_target, v_RC_source, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = l2-1-d;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
            }
            else            // r2->r1
            {
              std::pair<BoostSTREdge, bool> e_search = boost::add_edge(v_target, v_source, *p_graph_);
              if(e_search.second)
              {
                (*p_graph_)[e_search.first].p_A = c;
              }
              else
              {
                std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
              }
              e_search = boost::add_edge(v_RC_source, v_RC_target, *p_graph_);
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
        }
      }
    }
    while(std::getline(fin, line));

  }
  fin.clear();
  fin.close();
  return true;
}

void StrGraph::CondenseGraph()
{
  clock_t t;
  t = clock();
  int total_node = 0;
  for(auto it = boost::vertices(*p_graph_).first; it != boost::vertices(*p_graph_).second; ++it)
  {
    ++total_node;
    if(boost::in_degree(*it, *p_graph_) <= 0)
    { // orphant read is not seed
      if(boost::out_degree(*it, *p_graph_) > 0)   (*p_graph_)[*it].ff_seed = true;
    }
    else
    {
      if(boost::out_degree(*it, *p_graph_) != 1)  (*p_graph_)[*it].ff_seed = true;
      if(boost::in_degree(*it, *p_graph_) > 1)    (*p_graph_)[*it].ff_seed = true;
    }
  }
  t = clock() - t;
  std::cout << total_node << " vertices in total\n";
  std::cout << "find seed takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";


  std::list<cycle_s> to_add_cycle;
  int cnt_cycle_node = 0;
  double condeseSE = 0, condenseCycle = 0, findComponent = 0;
  BooleanVector connected = BooleanVector(2*p_order, false);
  for(auto it = boost::vertices(*p_graph_).first; it != boost::vertices(*p_graph_).second; ++it)
  {
    if(! p_crossed[ (*p_graph_)[*it].rid_ ]) // untraversed vertex, new component
    {
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
        if(!connected[rid])
        {
          connected[rid] = true;
          this_component.push_back(top_vertex);
          ++cnt_vertex;
          p_crossed[rid] = true;
          if(rid < p_order) p_crossed[rid+p_order] = true;
          else              p_crossed[rid-p_order] = true;
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
      if(ff_singleCycle)  // if this component is single loop
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
  std::cout << "find components takes " << findComponent << " s\n";
  std::cout << "consense SE takes " << condeseSE << " s\n";
  std::cout << "consense cycle takes " << condenseCycle << " s\n";
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
  std::cout << "modify graph takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
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
  clock_t t;
  t = clock();
  int cnt = 0;
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
  t = clock() - t;
  std::cout << "write "<< cnt << " edges " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  t = clock();
  cnt = 0;

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

  t = clock() - t;
  std::cout << "write "<< cnt << " orphants " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  fout.close();
}

*/


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
