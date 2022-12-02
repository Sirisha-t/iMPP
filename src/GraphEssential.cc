# include "GraphEssential.h"
#include<map>
#include<time.h>
using namespace std;

// prints the graph information
// parameters:
//    c: whether to print complete information including vertice and edge details; true for yes, false for no
void GraphEssential::PrintInfo(const bool c)  {
    std::cout << "Printing GraphEssential object info..." << std::endl;
    if(!initialized_)  {
        std::cout << "Graph is empty." << std::endl;
        return;
    }
    std::cout << "Num. of vertices:  " << boost::num_vertices(*graph_ptr_) << endl;
    std::cout << "Num. of edges:  " << boost::num_edges(*graph_ptr_) << endl;
    if(c) {
        // print complete information
        std::cout << "Printing all vertices information..." << std::endl;
        auto it_v = boost::vertices(*graph_ptr_).first;
        while(it_v != boost::vertices(*graph_ptr_).second) {
            std::cout << "===   begin of vertex:    " << std::endl;
            (*graph_ptr_)[*it_v].PrintInfo();
            std::cout << "===   end of vertex:    " << std::endl;
            ++ it_v;
        }
        std::cout << "Printing all edges information..." << std::endl;
        auto it_e = boost::edges(*graph_ptr_).first;
        while(it_e != boost::edges(*graph_ptr_).second) {
            std::cout << "===   begin of edge:    " << std::endl;
            (*graph_ptr_)[*it_e].PrintInfo();
            (*graph_ptr_)[boost::source(*it_e, *graph_ptr_)].PrintInfo();
            (*graph_ptr_)[boost::target(*it_e, *graph_ptr_)].PrintInfo();
            std::cout << "===   end of edge:  " << std::endl;
            ++ it_e;
        }
    }
    return;
}
// check whether the graph contains valid structure
bool GraphEssential::CheckGraphValidity(void)   {
    bool is_valid = true;
    bool check_passed = false;
    // checking overlap index bound
    check_passed = CheckOverlapIndex(graph_ptr_);
    is_valid = (is_valid && check_passed);
    // TODO: add other checks here
    return is_valid;
}

// Write the graph to a file with the ASQG format
void GraphEssential::WriteGraphASQG(const std::string & file)   {
    ofstream asqg_fh(file);
    if(!asqg_fh.is_open())   {
        cerr << "Error: GraphEssential::WriteGraphASQG: cannot open file to write graph. Abort." << endl;
        exit(1);
    }
    // writes the header info
    asqg_fh << "HT\tMANA Graph Module Output" << endl;
    // writes the vertex info
    auto it_v = boost::vertices(*graph_ptr_).first;
    while(it_v != boost::vertices(*graph_ptr_).second) {
        asqg_fh << "VT\t" << (*graph_ptr_)[*it_v].id_ << "\t" << (*graph_ptr_)[*it_v].str_ << "\tSS:i:0" << endl;
        ++ it_v;
    }
    // writes the edge info
    auto it_e = boost::edges(*graph_ptr_).first;
    while(it_e != boost::edges(*graph_ptr_).second) {
        BoostNodeType s = boost::source(*it_e, *graph_ptr_);
        BoostNodeType t = boost::target(*it_e, *graph_ptr_);
        asqg_fh << "ED\t" << (*graph_ptr_)[s].id_ << " " << (*graph_ptr_)[t].id_ << " "
            << (*graph_ptr_)[*it_e].ov_pos_[0] << " " << (*graph_ptr_)[*it_e].ov_pos_[1] << " " << (*graph_ptr_)[s].str_len_ << " "
            << (*graph_ptr_)[*it_e].ov_pos_[2] << " " << (*graph_ptr_)[*it_e].ov_pos_[3] << " " << (*graph_ptr_)[t].str_len_ << " "
            << (*graph_ptr_)[*it_e].is_rc_ << " " << (*graph_ptr_)[*it_e].num_mismatch_ << endl;
        ++ it_e;
    }
    asqg_fh.close();
    return;
}

void GraphEssential::RenameASQG(const std::string &file){
    ifstream asqg_fh (file);
    string newname = file;
    string outname1 = newname.substr(0,newname.find_last_of('.')) + ".rename.asqg";
    string outname2 = newname.substr(0,newname.find_last_of('.')) + ".namemap.txt";
    ofstream fout1(outname1.c_str());
    ofstream fout2(outname2.c_str());
    
    if(asqg_fh.is_open())    {
      boost::char_separator<char> sep(" \t");
      boost::tokenizer<boost::char_separator<char> >::iterator it;
      std::string line;
      vector<std::string>     p_header;
      unordered_map<std::string, int>  p_readID;
            
        while(getline(asqg_fh, line)) {
            if(line[0] == 'V')    {
              boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
              it = tokens.begin();
              fout1 << *it << "\t";
              it++;                           // cont[1], read name
              int id = p_readID.size();
              p_readID[*it] = id;
              fout1 << id << "\t";
              fout2 << *it << "\t"<<id<<"\n";
              it++;                           // cont[2], sequence
              fout1 << *it << "\n";
            }
            else if(line[0] == 'E'){
                boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
                it = tokens.begin();            // cont[0], ED
                fout1 << *it << "\t";
                it++;                           // cont[1], read1 name
                int source_id = p_readID[*it];
                fout1 << source_id << " ";
                it++;                           // cont[2], read2 name
                int target_id = p_readID[*it];
                fout1 << target_id;
                it++;
                while(it != tokens.end())
                {
                fout1 << " " << *it;
                it++;
                }
                fout1 << "\n";

            }
            else{
                fout1 << line << "\n";
            }
              
        }
        
    }
    asqg_fh.close();
    fout1.close();
    fout2.close();
}
  

void GraphEssential::LoadGraphASQG(const std::string &file)  {
    clock_t t= clock();
    IDType num_nodes = GetNumReadsASQG(file);
    t = clock() - t;
    std::cout<<"DEBUG: Number of reads in asqg file: "<<num_nodes<<"\n";
    cout<<"Time to GetNumReadsASQG :"<<((double)t)/CLOCKS_PER_SEC<<" s"<<endl;
    vector<bool> is_set(num_nodes, false);              // array indicating whether the orientation of a read is set
    vector<BoostNodeType> node_map(num_nodes);          // mapping from the node ID to the node decriptor
    BoostNodeType vt_node;
    unordered_map<IDType, BoostNodeType> vertex; // both
    //unordered_map< size_t, int> edgeCount_hash;
    //map< pair<int,int>, int> edgeCount_hash;

    size_t count=0;
    ifstream asqg_fh (file);
    if(asqg_fh.is_open())    {
        string line;
        t=clock();
        boost::char_separator<char> sep(" \t");
        boost::tokenizer<boost::char_separator<char> >::iterator it;
        while(getline(asqg_fh, line)) {

            //StringUtils::SplitByDelimiter(line, "\t ", vs);
            //assert(vs.size() >= 1);     //  fail assertion is likely due to sequence corruption
            if(line[0] == 'V')    {
              //vector<string> vs;
              boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
                // record the maximun read ID so far
                //assert(vs.size() >= 3);             // fail assertion sis likely due to sequence corruption
                it = tokens.begin();
                //IDType id = stoi(vs[1]);
                ++it;
                IDType id = stoi(*it);
                ++it;
                GraphNodeType *n = new GraphNodeType((*it).c_str());    // copying the sequence to the node; by default orientation = true (plus strand)
                n->id_ = id;                                            // DEBUG
                //cout << "DEBUG: raw sequence:   " << vs[2] << endl;
                //cout << "DEBUG: node sequence:   " << n.str_ << endl;
                //node_map[id] = boost::add_vertex(*n, *graph_ptr_);  // adding the node
                is_set[id] = true;

                /** ST test **/
                vt_node = boost::add_vertex(*n, *graph_ptr_);
                vertex[id] = vt_node;
            } else if(line[0] == 'E' )   {
              count++;
               //std::cout<<"ED :"<<count<<endl;
            //  if(count < 10000000)
            //  {

               boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
               it = tokens.begin();
               ++it;
                t = clock();
                IDType a = stoi(*it);     // the source ID
                ++it;
                IDType b = stoi(*it);     // the target ID
                // checking validity of the edge
                if(a == b)  {               // a self-loop
                    cerr << "Warning:   GraphEssential::ReadASQG: A self-loop detected, edge ignored." << endl;
                    continue;
                }   else if (!is_set[a] || !is_set[b])  {   // one of the nodes does not exist in the graph
                    cerr << "Warning:   GraphEssential::ReadASQG: Trying to add edge to nodes that do not exist, edge ignored." << endl;
                }
                //cout << "DEBUG: two vertices:   " << a << " " << b << endl;
                // determine which node is the source and which is the target
                IDType s, t;
                bool src_ori;   // orientation of the source read
                bool olp_rev;   // whether the overlap positions should be reversed
                ++it;
                int s1 = std::stoi(*it);
                ++it;
                int e1 = std::stoi(*it);
                ++it;
                int l1 = std::stoi(*it);
                ++it;
                int s2 = std::stoi(*it);
                ++it;
                int e2 = std::stoi(*it);
                ++it;
                int l2 = std::stoi(*it);
                ++it;
                int rc = std::stoi(*it) ;
                ++it;
                int mm = std::stoi(*it);
                if(s1 == 0)    {
                    // the first read has its prefix overlapped
                    t = a; s = b;
                    src_ori = (bool) (rc) ? false : true;    // a is the target; if the overlap is reverse complementary then the source is negative strand; otherwise positive
                    olp_rev = true;
                }   else   {
                    // the first read has its suffix overlapped
                    s = a; t = b;
                    src_ori = true;     // always be true, because a is always positive strand
                    olp_rev = false;
                }
                //std::cout<<a<<" "<<" "<<b<<" "<<s1<<" "<<s2<<" "<<e1<<" "<<e2<<" "<<rc<<"\n";
                /** ST test **/
                //vt_source = boost::add_vertex(*n, *graph_ptr_);
                //vertex[s] = vt_source;
                //vt_target = boost::add_vertex(*n, *graph_ptr_);
                //vertex[t] = vt_target;


                // Check if an edge already exists
                //pair<BoostEdgeType, bool> e_check = boost::edge(node_map[s], node_map[t], *graph_ptr_);  // note that the target is the first parameter and the source is the second
                //if(!e_check.second)    {
                    //cout << "DEBUG: edge does not exist." << endl;
                    pair<BoostEdgeType, bool> e_add = boost::add_edge(vertex[s], vertex[t], *graph_ptr_);     // note that the target is the first parameter and the source is the second
                    if(!e_add.second)    {
                        cerr << "Warning:   GraphEssential::ReadASQG: Failed to add edge, edge ignored." << endl;
                    }   else    {
                        // if edge is successfully added, incorate related information
                        if((*graph_ptr_)[e_add.first].GetOverlap() < (e1 - s1 + 1))
                          (*graph_ptr_)[e_add.first].SetOverlap(e1 - s1 + 1);
                        (*graph_ptr_)[e_add.first].SetIsRevComplement((bool) rc);  // "1" in the ASQG file indicates reverse complement
                        (*graph_ptr_)[e_add.first].SetSrcOrientation(src_ori);
                        (*graph_ptr_)[e_add.first].SetNumMismatches(mm);
                        if(!olp_rev)    {
                            (*graph_ptr_)[e_add.first].SetOverlapPosition(s1, e1, s2, e2);
                        }   else    {
                            (*graph_ptr_)[e_add.first].SetOverlapPosition(s2, e2, s1, e1);
                        }
                        //cout << "DEBUG: source length:  " << (*graph_ptr_)[node_map[s]].str_len_ << endl;
                        //cout << "DEBUG: target length:  " << (*graph_ptr_)[node_map[t]].str_len_ << endl;
                    }
                //}
                   /*else    {
                    // check overlap, use the larger overlap between sequences
                    //cout << "DEBUG: duplicate edge overlap info:    " << (*graph_ptr_)[e_check.first].GetOverlap() << " " << stoi(vs[4]) - stoi(vs[3]) + 1 << endl;
                    //cout << "DEBUG: duplicate edge rev_comp info:    " << (*graph_ptr_)[e_check.first].IsRevComplement() << " " << stoi(vs[9]) << endl;
                    if((*graph_ptr_)[e_check.first].GetOverlap() < stoi(ed[4]) - stoi(ed[3]) + 1)    {
                        (*graph_ptr_)[e_check.first].SetOverlap(stoi(ed[4]) - stoi(ed[3]) + 1);
                        (*graph_ptr_)[e_check.first].SetIsRevComplement((bool) stoi(ed[9]));  // "1" in the ASQG file indicates reverse complement
                        (*graph_ptr_)[e_check.first].SetSrcOrientation(src_ori);
                        (*graph_ptr_)[e_check.first].SetNumMismatches(stoi(ed[10]));
                        if(!olp_rev)    {
                            (*graph_ptr_)[e_check.first].SetOverlapPosition(stoi(ed[3]), stoi(ed[4]), stoi(ed[6]), stoi(ed[7]));
                        }   else    {
                            (*graph_ptr_)[e_check.first].SetOverlapPosition(stoi(ed[6]), stoi(ed[7]), stoi(ed[3]), stoi(ed[4]));
                        }
                    }
                }*/
            //}else{ break;}

          //}while(getline(asqg_fh, line));
          //else{
          //  break;
          //}
        }
      }
        t = clock() - t;
        cout<<"Time to read ED (vector)(ST): "<<((double)t)/CLOCKS_PER_SEC <<"s"<< endl;
        cout<<"Edges read : "<<count<<endl;

          asqg_fh.close();
        }
        return;
      }


IDType GraphEssential::GetNumReadsASQG(const std::string &file) {
    IDType max_ID = 0;
    ifstream asqg_fh (file);
    if(asqg_fh.is_open())    {
      boost::char_separator<char> sep(" \t");
      boost::tokenizer<boost::char_separator<char> >::iterator it;
        string line;
        while(getline(asqg_fh, line)) {
            vector<string> vs;
            //StringUtils::SplitByDelimiter(line, "\t ", vs);
            //assert(vs.size() >= 1);     //  fail assertion is likely due to sequence corruption
            if(line[0] == 'V')    {
              boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
              it = tokens.begin();
              ++it;
                // record the maximun read ID so far
                //assert(vs.size() >= 3);     // fail assertion is likely due to sequence corruption
                max_ID = max_ID < stoi(*it) ? stoi(*it) : max_ID;
            }   else if(line[0] == 'E')   {
                break;
            }
        }
        asqg_fh.close();
    }
    return (max_ID + 1);   // we need to add 1 to convert the ID into the size
}

// check whether the overlap indexes are out of bound
bool GraphEssential::CheckOverlapIndex(AssemblyGraphType *g)    {
    auto it_e = boost::edges(*g).first;
    while(it_e != boost::edges(*g).second) {
        BoostNodeType s = boost::source(*it_e, *g);
        BoostNodeType t = boost::target(*it_e, *g);
        // check the source positions
        if((*g)[*it_e].ov_pos_[0] > (*g)[s].str_len_ - 1 || (*g)[*it_e].ov_pos_[1] > (*g)[s].str_len_ - 1)  {
            cerr << "GraphEssential::CheckOverlapIndex: source overlap index out of bound. "
                    << (*g)[s].id_ << "  " << (*g)[t].id_ << "  "
                    << (*g)[*it_e].ov_pos_[0] << "    " << (*g)[*it_e].ov_pos_[1] << "  "
                    << (*g)[s].str_len_ << endl;
            return false;
        }

        // check the target positions
        if((*g)[*it_e].ov_pos_[2] > (*g)[t].str_len_ - 1 || (*g)[*it_e].ov_pos_[3] > (*g)[t].str_len_ - 1)  {
            cerr << "GraphEssential::CheckOverlapIndex: target overlap index out of bound. "
                    << (*g)[s].id_ << "  " << (*g)[t].id_ << "  "
                    << (*g)[*it_e].ov_pos_[2] << "    " << (*g)[*it_e].ov_pos_[3] << "  "
                    << (*g)[t].str_len_ << endl;
            return false;
        }
        ++ it_e;
    }
    return true;
}
