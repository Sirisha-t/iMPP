#include "SG.h"
#include "unionFind.hpp"


void StrGraph::writeMergedGraph(std::string & filename)
{
  /** timer **/
  clock_t t;
  t = clock();
  /***********/
  std::cout << "start write..\n";
  std::string outname = filename.substr(0, filename.find_last_of('.')) + ".merged.fq";
  std::ofstream fout(outname.c_str());
  for(int node_i = 0; node_i < p_seq.size(); node_i++)
  {
    fout << ">" << node_i << ",";
    for(int r = 0; r < p_read_on_node[node_i].size(); r++)
    {
      fout << p_read_on_node[node_i][r] << ',' << p_pos_on_node[node_i][r] << ',';
      if(p_read_on_node[node_i][r] < 0)        fout << "-1";
      else                                     fout << p_read_length[p_read_on_node[node_i][r]];
      if(r+1 < p_read_on_node[node_i].size())  fout << ',';
    }

    BoostSTRVertex v = vertex[node_i];
    auto it_v = boost::adjacent_vertices(v, *p_graph_).first;
    if(it_v != boost::adjacent_vertices(v, *p_graph_).second)
    {
      fout << ':';
      for(; it_v != boost::adjacent_vertices(v, *p_graph_).second; it_v++)
      {
        fout << (*p_graph_)[*it_v].rid_;
        it_v++;
        if(it_v != boost::adjacent_vertices(v, *p_graph_).second) fout << ',';
        it_v--;
      }
    }
    fout << "\n" << p_seq[node_i] << "\n";
  }
  /** timer **/
  t = clock() - t;
  std::cout << "write graph takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  /***********/
  fout.close();
}

void StrGraph::mergeSG(std::string& sam_name, std::string& contig_name)
{
  /*********/
  std::ifstream fin;
  /*********/
  clock_t t;
  /*********/
  /* load contigs */
  /*********/
  t = clock();
  /*********/
  fin.open(contig_name.c_str());
  StringVector contigs;
  int count = 0;
  String2Integer contigs_id;
  if (!fin.is_open())	    // fail to open
  {
    std::cout << "Error: " << contig_name << " doesnot exist!\n";
  }
  else				// succeed to open
  {
    std::string line;
    std::string name;
    while(std::getline(fin, line))
    {
      if(line[0]=='>')
      {     
        name = line.substr(1);
        contigs_id[name] = contigs_id.size();    
       
      }
      else
      {
        if(contigs.size() < contigs_id.size())  contigs.push_back("");
        contigs[contigs_id.size()-1] += line;
       
      }
    }
  }
  fin.clear();
  fin.close();
  /*********/
  t = clock() - t;
  std::cout << "Loading contigs takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  /*********/
  /* load terminals */
  /*********/
  t = clock();
  /*********/
  std::vector<BooleanVector> taken_2D;
  for(auto x : contigs) taken_2D.push_back(BooleanVector(x.size(), false));
  std::vector<BooleanVector> coverage_2D;
  for(auto x : contigs) coverage_2D.push_back(BooleanVector(x.size(), false));
  std::vector<std::map<int, int>> f_terminals = std::vector<std::map<int, int>>(contigs.size());
  std::vector<std::map<int, int>> f_match = std::vector<std::map<int, int>>(contigs.size());
  fin.open(sam_name.c_str());	// open file
  if (!fin.is_open())	    // fail to open
  {
    std::cout << "Error: " << sam_name << " does not exist!\n";
  }
  else				// succeed to open
  {
    std::string line;
    boost::char_separator<char> sep("\t");
    boost::tokenizer<boost::char_separator<char> >::iterator it;
    while(std::getline(fin, line))
    {
      if(line[0] != '@')
      {
        boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
        it = tokens.begin();  // edge name
        int node_id = std::stoi((*it));
        std::string header = (*it);
        bool ff_orphant_read = thisIsOrphat( header );
        ++it;                 // flag
        int flag = std::stoi((*it));
        ++it;                 // reference name
        std::string cont_name = *it;
        int ref_id = contigs_id[(*it)];
        ++it;                 // pos
        std::string pos_str = (*it);
        ++it;                 // MAPQ
        ++it;                 // CIAGR
        std::string CIAGR = (*it);
        BoostSTRVertex node = vertex[node_id];
        //std::cout<<node_id<<"\t"<<header<<"\t"<<ff_orphant_read<<"\t"<<flag<<"\t"<<cont_name<<"\t"<<ref_id<<"\t"<<pos_str<<"\t"<<CIAGR<<"\t"<<node<<"\n";

        
        unsigned int mapping_quality = goodMapping(CIAGR, flag);
        if(ff_orphant_read)
        {
          if((mapping_quality&0xF0000000) == 0x00000000)
          {
            int read_length = p_seq[node_id].size();
            int pos = std::stoi(pos_str);
            for(int i = 0; i < read_length; ++i)
            {
              coverage_2D[ref_id][pos+i-1] = true;
              taken_2D[ref_id][pos+i-1] = true;
            }
          }
        }
        else
        {
          int pos = std::stoi(pos_str);
          int match_length = (mapping_quality&0x0FFFFFFF);
          for(int i = 0; i < match_length; ++i) taken_2D[ref_id][pos+i-1] = true;
          if(((mapping_quality&0xE0000000) == 0x00000000) && ((boost::in_degree(node, *p_graph_)*boost::out_degree(node, *p_graph_)) == 0))
          {
            if((flag&16) == 0)
            {
              f_terminals[ref_id][pos] = node_id;
              f_match[ref_id][node_id] = match_length;
            }
          }
        }

      }
    }
  }
  fin.clear();
  fin.close();
  /*********/
  t = clock() - t;
  std::cout << "load terminals takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  /*********/

  /* building DisjointSet */
  /*********/
  t = clock();
  /*********/
  DisjointSet* ds = new DisjointSet();
  for(auto it_e = boost::edges(*p_graph_).first; it_e != boost::edges(*p_graph_).second; ++it_e)
  {
    BoostSTRVertex head = boost::source(*it_e, *p_graph_);
    BoostSTRVertex tail = boost::target(*it_e, *p_graph_);
    ds->unionNode((*p_graph_)[head].rid_, (*p_graph_)[tail].rid_);
  }
  /*********/
  t = clock() - t;
  std::cout << "construct DisjointSet takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  /*********/

  /*********/
  t = clock();
  /*********/
  for(int ref_id = 0; ref_id < f_terminals.size(); ref_id++)
  {
    auto this_terminal = f_terminals[ref_id];
    auto it1 = this_terminal.begin();
    auto it2 = this_terminal.begin();
    std::string seq = contigs[ref_id];

    // take orphant seq
    int p = 0, q = 0;
    while(p < seq.size())
    {
      while((p < seq.size()) && taken_2D[ref_id][p]) p++;
      if(p - q > 200)
      {
        std::string new_seq = seq.substr(q, p-q);
        IntegerVector read_on_this_node, loc_on_this_node;
        read_on_this_node.push_back(1);
        loc_on_this_node.push_back(1);

        int new_node_id = p_seq.size();
        p_seq.push_back(new_seq);
        p_read_on_node.push_back(read_on_this_node);
        p_pos_on_node.push_back(loc_on_this_node);

        BoostSTRVertex new_node;
        if(vertex.count(new_node_id) == 0)
        {
          STRVertexType node(new_node_id);
          new_node = boost::add_vertex(node, *p_graph_);
          vertex[new_node_id] = new_node;
        }
        else
        {
          new_node = vertex[new_node_id];
        }
      }
      p++;
      q = p;
    }

    // merge
    if(it2 != this_terminal.end())
    {
      // before first
      int pos1 = it1->first;
      int node1 = it1->second;
      int s = pos1-2;
      while((s>=0) && coverage_2D[ref_id][s]) s--;
      if(s != (pos1-2))
      {
        int first_read_id = p_read_on_node[node1][0];
        int len = pos1 + p_read_length[first_read_id] - s - 2;
        //std::cout<<seq<<"\t"<<s+1<<"\t"<<len<<"\n";
        std::string new_seq = seq.substr(s+1, len);
	//std::cout<<new_seq<<"\n";
        IntegerVector read_on_this_node, loc_on_this_node;
        read_on_this_node.push_back(first_read_id);
        loc_on_this_node.push_back(pos1 - s + 1);

        int new_node_id = p_seq.size();
        p_seq.push_back(new_seq);
        p_read_on_node.push_back(read_on_this_node);
        p_pos_on_node.push_back(loc_on_this_node);

        BoostSTRVertex new_node;
        if(vertex.count(new_node_id) == 0)
        {
          STRVertexType node(new_node_id);
          new_node = boost::add_vertex(node, *p_graph_);
          vertex[new_node_id] = new_node;
        }
        else
        {
          new_node = vertex[new_node_id];
        }

        BoostSTRVertex v_target = vertex[node1];
        std::pair<BoostSTREdge, bool> e_search2 = boost::add_edge(new_node, v_target, *p_graph_);
        if(!e_search2.second)          std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
      }

      it2++;
      for(; it2 != this_terminal.end(); it1++, it2++)
      {
        int pos1 = it1->first;
        int node1 = it1->second;
        int pos2 = it2->first;
        int node2 = it2->second;
        int end = pos1 + f_match[ref_id][node1] - 2;              // coordinate of last nt of node1 on contig
        if(ds->find(node1) != ds->find(node2)) // node1 and node2 should not be connected!!!
        {
          int last_read_id = *(p_read_on_node[node1].end()-1);
          int first_read_id = p_read_on_node[node2][0];
          int l1 = p_read_length[last_read_id];
          int l2 = p_read_length[first_read_id];
          int y = end + 1 - l1;
          int z = pos2 + l2 -2;
          int len = z - y + 1;

          if(len == l1)
          {
            BoostSTRVertex v_source = vertex[node1];
            BoostSTRVertex v_target = vertex[node2];
            std::pair<BoostSTREdge, bool> e_search1 = boost::add_edge(v_source, v_target, *p_graph_);
            if(!e_search1.second)          std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
          }
          else
          {
            std::string new_seq;
            IntegerVector read_on_this_node, loc_on_this_node;
            if(len > l1)
            {
              new_seq = seq.substr(y, len);
              read_on_this_node.push_back(last_read_id);
              read_on_this_node.push_back(first_read_id);
              loc_on_this_node.push_back(1);
              loc_on_this_node.push_back(pos2 - y);
            }
            else
            {
              int cut = end - z + 1;
              new_seq = "N";
              read_on_this_node.push_back(-1);
              loc_on_this_node.push_back(-cut);
            }

            // add node
            int new_node_id = p_seq.size();
            p_seq.push_back(new_seq);
            p_read_on_node.push_back(read_on_this_node);
            p_pos_on_node.push_back(loc_on_this_node);

            BoostSTRVertex new_node;
            if(vertex.count(new_node_id) == 0)
            {
              STRVertexType node(new_node_id);
              new_node = boost::add_vertex(node, *p_graph_);
              vertex[new_node_id] = new_node;
            }
            else
            {
              new_node = vertex[new_node_id];
            }

            BoostSTRVertex v_source = vertex[node1];
            std::pair<BoostSTREdge, bool> e_search1 = boost::add_edge(v_source, new_node, *p_graph_);
            if(!e_search1.second)          std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";

            BoostSTRVertex v_target = vertex[node2];
            std::pair<BoostSTREdge, bool> e_search2 = boost::add_edge(new_node, v_target, *p_graph_);
            if(!e_search2.second)          std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
          }
        }
      }
      // after last
      pos1 = it1->first;
      node1 = it1->second;
      int e = pos1 + f_match[ref_id][node1]-1;
      while((e < coverage_2D[ref_id].size()) && coverage_2D[ref_id][e]) e++;
      if(e != (pos1 + f_match[ref_id][node1]-1))
      {
        int last_read_id = p_read_on_node[node1][0];
        int s = pos1 + f_match[ref_id][node1] - p_read_length[last_read_id] - 1;
        int len = e - s;
        std::string new_seq = seq.substr(s, len);
        IntegerVector read_on_this_node, loc_on_this_node;
        read_on_this_node.push_back(last_read_id);
        loc_on_this_node.push_back(1);

        int new_node_id = p_seq.size();
        p_seq.push_back(new_seq);
        p_read_on_node.push_back(read_on_this_node);
        p_pos_on_node.push_back(loc_on_this_node);

        BoostSTRVertex new_node;
        if(vertex.count(new_node_id) == 0)
        {
          STRVertexType node(new_node_id);
          new_node = boost::add_vertex(node, *p_graph_);
          vertex[new_node_id] = new_node;
        }
        else
        {
          new_node = vertex[new_node_id];
        }

        BoostSTRVertex v_source = vertex[node1];
        std::pair<BoostSTREdge, bool> e_search1 = boost::add_edge(v_source, new_node, *p_graph_);
        if(!e_search1.second)          std::cout << "Error:: StringGraph::LoadGraph: Failed to add edges between vertices!\n";

      }
    }
  }
  /*********/
  t = clock() - t;
  std::cout << "merge graph takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  /*********/

  delete ds;
}


bool thisIsOrphat(std::string& header)
{
  boost::char_separator<char> sep(",");
  boost::tokenizer<boost::char_separator<char> > tokens(header, sep);
  boost::tokenizer<boost::char_separator<char> >::iterator it;
  it = tokens.begin();
  ++it; // discard read_id
  int i;
  for(i = 0 ; it != tokens.end(); ++it, ++i) if(i > 3) return false;
  return (i == 3);
}

unsigned int goodMapping(std::string& CIAGR, int flag)
{
  if((flag&4) == 4)                 return 0x80000000;  // unmap
  unsigned int match = 0;
  std::string temp = "";
  std::string symbols = "";
  for(int i = 0; i < CIAGR.size(); i++)
  {
    if(CIAGR[i] >= '0' && CIAGR[i] <= '9')    temp += CIAGR[i];
    else
    {
        symbols += CIAGR[i];
        if(CIAGR[i] == 'M') match += std::stoi(temp);
        temp = "";
    }
  }
  if(match < 100)                   return 0x20000000;  // mapping length too short
  if(symbols.size() > 2)            return 0x20000000;  // too fragmented
  if (symbols.size() == 2)
  {
    if(((flag&16) == 16) && symbols[1] == 'M')
                                    return 0x40000000;  // tail unmap
    if(((flag&16) == 0)  && symbols[1] != 'M')
                                    return 0x40000000;  // tail unmap
    return (0x10000000+match);
  }
  return match;
}


bool StrGraph::loadAnchor(
  char* pEdgeSet,
  int MAX_EXTEND_LENGTH,
  anchors& anchor_set)
{
  std::ifstream fin;
  fin.open(pEdgeSet);
  if (!fin.is_open())
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

    int length = p_seq[node].size();
    //std::cout<<"Length of node is: "<<length<<"\n";
    int left, right;
    //new_anchor.push_back(node);
    if(MAX_EXTEND_LENGTH > 0)
    {
        left = (MAX_EXTEND_LENGTH - length)/2;
        right = (MAX_EXTEND_LENGTH - length)/2;
    }
    new_anchor.path = IntegerList{node};
    new_anchor.M_upstream = left;
    new_anchor.M_downstream = right;


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

bool StrGraph::extendGraph(
  char* fgsoutname,
  int MAX_EXTEND_LENGTH,
  char* path_name,
  bool db_flag)
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
  p_crossed_h = BooleanVector(p_seq.size(), false);
  p_crossed_t = BooleanVector(p_seq.size(), false);
  int anc_count=0;
  size_t total_path = 0;

  //std::cout<<"Anchors:"<<anchor_set.size()<<"\n";
  for(auto it = anchor_set.begin(); it != anchor_set.end(); ++it)
  {
    int anc_path=0;
    anc_count++;
    //std::cout<<"Anchor to be traversed: "<<anc_count<<"\n";
    DirectedDFS((*it),prededge_set, MAX_EXTEND_LENGTH, path_name, anc_path, db_flag);
    size_t total_path = total_path + anc_path;
  }

  /** timer **/
  t = clock() - t;
  std::cout << "Extending anchor takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  /***********/
  return true;
}

void StrGraph::DirectedDFS(struct anchor& anchor, predEdges& prededge_set, int MAX_EXTEND_LENGTH, char* path_name, int& pathnum, bool &ff_spadesGraph)
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
            //p_crossed[new_head] = true;
            int last_read_length;
	    if(!p_crossed_h[new_head]){
              p_crossed_h[new_head] = true;
            if(prededge_set.count(new_head) == 0) //unpredicted adjacent edge
            {
              if(!ff_spadesGraph){
                IntegerVector reads = p_read_on_node[new_head];
                int last_read_id = *(reads.end()-1);
                last_read_length = p_read_length[last_read_id];
              }
              else{
                last_read_length = kmer_length;
              }

              // update anchor stack
              struct anchor new_anchor;
              new_anchor.path = path;
              new_anchor.M_upstream = top_anchor.M_upstream - p_seq[new_head].size() + last_read_length; // +last_read_length;
              new_anchor.M_downstream = top_anchor.M_downstream;
              //std::cout<<"head :"<<new_head<<"\t"<<"upstream :"<<new_anchor.M_upstream<<"\n";

              checking_stacks.push(new_anchor);
              path.pop_front();
            }
            else //predicted edge
            {
		    //std::cout<<"Predicted adj edge(up) : "<<new_head<<"\n";
                    top_anchor.path = path;
                    half_done_stacks.push(top_anchor);
                    path.pop_front();

            }
	}
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
        if(it_v == boost::adjacent_vertices(last_node, *p_graph_).second)  // no tail
        {
                 pathnum++;
                 if(ff_spadesGraph) savePath(path_name, top_anchor.path);
                 else               savePath(path_name, top_anchor);
        }
        else
        {
          for(; it_v != boost::adjacent_vertices(last_node, *p_graph_).second; ++it_v)
          {
            int new_tail = (*p_graph_)[*it_v].rid_;
            if(!formCycle(path, new_tail))  // no cycle allowed
            {
              path.push_back(new_tail);
              //p_crossed[new_tail] = true;
              int first_read_length;
              if( !p_crossed_t[new_tail]){
                p_crossed_t[new_tail] = true;
	      if(prededge_set.count(new_tail) == 0)
              {
                if(!ff_spadesGraph)
                {
                  IntegerVector reads = p_read_on_node[new_tail];
                  int first_read_id = reads[0];
                  first_read_length = p_read_length[first_read_id];
                }
                else
                {
                  first_read_length = kmer_length;
                }
                // update anchor stack
                struct anchor new_anchor;
                new_anchor.path = path;
                new_anchor.M_upstream = top_anchor.M_upstream;
                new_anchor.M_downstream = top_anchor.M_downstream - p_seq[new_tail].size() + first_read_length;

                half_done_stacks.push(new_anchor);
                path.pop_back();
              }
              else{
		    //std::cout<<"Pred adj edge(down) : "<<new_tail<<"\n";
                    top_anchor.path = path;
                    if(ff_spadesGraph)    savePath(path_name, top_anchor.path);
                    else                  savePath(path_name, top_anchor);
                    path.pop_back();
              }
	    }

            }
        }
      }
    }
    else{
          pathnum++;

          if(ff_spadesGraph)    savePath(path_name, top_anchor.path);
          else                  savePath(path_name, top_anchor);
    }
  }

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
    if(it != path.begin()) 
	{ 
	  if(first_read < 0) //merged node
		seq = seq.substr(0, seq.size() + p_pos_on_node[edge][0]);
    	  else	seq += full.substr(p_read_length[first_read]);
	}
    else
	seq += full;
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

// Reading string graph file input
bool StrGraph::readSGFile(std::string & filename)
{
  /** timer **/
  clock_t t;
  t = clock();
  /***********/
  std::ifstream fin;
  // fin.open(config.IDIR_GRAPH.c_str());	// open file
  fin.open(filename.c_str());	// open file
  if (!fin.is_open())	    // fail to open
  {
    std::cout << "Error: " << filename << " doesnot exist!\n";
    return false;
  }
  else				// succeed to open
  {
    std::cout<<"Reading String Graph file .. \n";
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

//Reading de Bruijn graph file Input
bool StrGraph::readSPAdesFile(std::string & filename, int km)
{
  /** timer **/
  clock_t t;
  t = clock();
  /***********/
  std::ifstream fin;
  // fin.open(config.IDIR_GRAPH.c_str());       // open file
  fin.open(filename.c_str());   // open fileq
  if (!fin.is_open())       // fail to open
  {
    std::cout << "Error: " << filename << " does not exist!\n";
    return false;
  }
  else                          // succeed to open
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
          int node_id = std::stoi(header);
	  p_all_nodes.push_back(node_id);

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
          int node_id = std::stoi(header);
          p_all_nodes.push_back(node_id);

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
  fin.clear();
  fin.close();
  p_traversed = BooleanVector(p_seq.size(), false);
  ff_spadesGraph = true;
  kmer_length = km;
  /** timer **/
  t = clock() - t;
  std::cout << "Loading de Bruijn graph takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
  /***********/
  return true;
}


void StrGraph::extendSPAdes(char* path_name, int MAXLEN)
{
  /** timer **/
  clock_t t;
  t = clock();
  /***********/
  std::ofstream path_out;
  path_out.open(path_name);
  path_out.clear();
  path_out.close();
  p_node_id = 0;
  // find all source vertex
  std::stack<IntegerList> to_visit;
  for(auto it = boost::vertices(*p_graph_).first; it != boost::vertices(*p_graph_).second; ++it)
  {
    if(boost::in_degree(*it, *p_graph_) == 0)
    {
      IntegerList init_path;
      int current_id = (*p_graph_)[*it].rid_;
      init_path.push_back(current_id);
      p_traversed[current_id] = true;
      if(boost::out_degree(*it, *p_graph_) > 0)
      {
        for(auto it_v = boost::adjacent_vertices(*it, *p_graph_).first; it_v != boost::adjacent_vertices(*it, *p_graph_).second; ++it_v)
        {
          int t_id = (*p_graph_)[*it_v].rid_;
          p_traversed[t_id] = true;
          init_path.push_back(t_id);
          to_visit.push(init_path);
          init_path.pop_back();
        }
      }
      else savePath(path_name, init_path);
    }
  }
  while(!to_visit.empty())
  {
    IntegerList top_path = to_visit.top();
    to_visit.pop();
    int last_node = top_path.back();
    BoostSTRVertex last_vertex = vertex[last_node];

   if(getPathLength(top_path) > MAXLEN)
   {
     savePath(path_name, top_path);
     if(boost::out_degree(last_vertex, *p_graph_) > 0)
     {
       IntegerList init_path{last_node};
       for(auto it_v = boost::adjacent_vertices(last_vertex, *p_graph_).first; it_v != boost::adjacent_vertices(last_vertex, *p_graph_).second; ++it_v)
       {
         int t_id = (*p_graph_)[*it_v].rid_;
         if(!p_traversed[t_id])
         {
           p_traversed[t_id] = true;
           init_path.push_back(t_id);
           to_visit.push(init_path);
           init_path.pop_back();
         }
       }
     }
   }
   else
   {
     if(boost::out_degree(last_vertex, *p_graph_) > 0)
     {
       for(auto it_v = boost::adjacent_vertices(last_vertex, *p_graph_).first; it_v != boost::adjacent_vertices(last_vertex, *p_graph_).second; ++it_v)
       {
         int t_id = (*p_graph_)[*it_v].rid_;
         if(!p_traversed[t_id])
         {
           p_traversed[t_id] = true;
           top_path.push_back(t_id);
           to_visit.push(top_path);
           top_path.pop_back();
         }
       }
     }
     else
     {
       savePath(path_name, top_path);
     }
   }
 }
   /** timer **/
   t = clock() - t;
   std::cout << "Extending SPAdes edges takes " << ((double)t)/CLOCKS_PER_SEC << " s\n";
   /***********/

 }

 int StrGraph::getPathLength(const IntegerList& path)
{
  int ans = kmer_length;
  for(auto it = path.begin(); it != path.end(); ++it)
  {
    ans += p_seq[*it].size() - kmer_length;
  }
  return ans;
}

void StrGraph::savePath(char* path_name, const IntegerList& path)
{
  ++p_node_id;
  std::ofstream path_out;
  path_out.open(path_name,std::ios::in|std::ios::out|std::ios::app);
  std::string seq;
  path_out << ">" << p_node_id;
  for(auto it = path.begin(); it != path.end(); ++it)
  {
    path_out << "," << *it;
    std::string full = p_seq[*it];
    if(it != path.begin()) seq += full.substr(kmer_length);
    else                   seq += full;
  }
  path_out << "\n" << seq << "\n";
  path_out.close();
}


StrGraph::StrGraph()
{
  p_graph_ = new BoostSTRGraph();
}

StrGraph::~StrGraph()
{
  delete p_graph_;
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
