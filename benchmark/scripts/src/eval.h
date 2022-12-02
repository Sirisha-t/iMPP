#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <time.h>
#include <regex>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <set>
#include <map>
#include <unordered_map>
#include <boost/tokenizer.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <time.h>
#include <boost/bimap.hpp>

/************************************************/
typedef std::vector<std::string> StringArray;
typedef boost::bimap<std::string, int> StringIntegerMap;
typedef std::unordered_map<int, float> IntegerFloatMap;
typedef std::map<std::string, int> StringMap;
typedef std::vector< std::vector<int> > Integer2D;
typedef std::vector<int> IntegerArray;
typedef std::set<int> IntegerSet;
typedef std::map<int, int> IntegerMap;
typedef std::map<int, IntegerSet> MapSet;

typedef std::pair<std::size_t, std::size_t> KeyPair;
typedef std::vector<KeyPair> VectorPair;
typedef std::map<int, VectorPair> PairMap;
/***********************************************/
double threshold = 0.6;
double rthresh = 0.6;
int freq = 40; //for ROC freq points
StringIntegerMap readNames; //bimap to store read names <-> number
clock_t begin, end, elapsed;
/*********************/
//File name variables
//bwa files
std::string bwa_ref = "";
std::string bwa_path = "";
std::string bwa_edge = "";
std::string bwa_sgacg = "";
std::string bwa_spadescg = "";
// FGS files
std::string fgs_ref = "";
std::string fgs_read = "";
std::string fgs_path = "";
std::string fgs_edge = "";
std::string fgs_sgacg = "";
std::string fgs_spadescg = "";

std::string reads_fastq="";
//6f
std::string sixfr_reads = "";
//Output eval files
std::string eval_read = "";
std::string eval_path = "";
std::string eval_edge = "";
std::string eval_sgacg = "";
std::string eval_spadescg = "";

//roc output files
std::string roc_read = "";
std::string roc_path = "";
std::string roc_edge = "";
std::string roc_sgacg = "";
std::string roc_spadescg = "";

/************************/
//Storing bwa ref
StringMap refName;
//PairMap read_RefQuery;
Integer2D read_RefQuery;
Integer2D pos_of_readRef;
//Integer2D spos_of_readRef;
//Integer2D epos_of_readRef;
IntegerMap lenRef;
IntegerArray refMapped;

//storing bwa edge
StringMap edgeName;
Integer2D read_EdgeQuery;
Integer2D pos_of_readEdge;
IntegerMap lenEdge;


// storing bwa path
StringMap pathName;
Integer2D read_PathQuery;
Integer2D pos_of_readPath;
IntegerMap lenPath;


//storing bwa sga contig
StringMap SGAcontigName;
Integer2D read_SGAContigQuery;
Integer2D pos_of_SGAreadContig;
IntegerMap lenSGACG;


//storing bwa spades contig
StringMap SPAcontigName;
Integer2D read_SPAContigQuery;
Integer2D pos_of_SPAreadContig;
IntegerMap lenSPACG;


/**********************/
Integer2D fgs_readsPerRef;
IntegerArray fgs_allReadsRef;


Integer2D fgs_readsPerRead;
Integer2D fgs_allReads;
IntegerArray fgs_allReadsArray;

Integer2D fgs_readsPerEdge;
Integer2D fgs_allReadsEdge;
IntegerArray fgs_allReadsEdgeArray;
MapSet predReadsEdge;


Integer2D fgs_readsPerPath;
Integer2D fgs_allReadsPath;
IntegerArray fgs_allReadsPathArray;
MapSet predReadsPath;

Integer2D fgs_readsPerSGACG;
Integer2D fgs_allReadsSGACG;
IntegerArray fgs_allReadsSGAArray;
MapSet predReadsSGA;

Integer2D fgs_readsPerSPACG;
Integer2D fgs_allReadsSPACG;
IntegerArray fgs_allReadsSPAArray;
MapSet predReadsSPA;

IntegerArray sixf_reads;

/**********************/
//storing tp, fp and fn
IntegerArray read_tp;
IntegerArray read_fp;
IntegerArray read_fn;

IntegerArray edge_tp;
IntegerArray edge_fp;
IntegerArray edge_fn;

IntegerArray path_tp;
IntegerArray path_fp;
IntegerArray path_fn;

IntegerArray sgacg_tp;
IntegerArray sgacg_fp;
IntegerArray sgacg_fn;

IntegerArray spacg_tp;
IntegerArray spacg_fp;
IntegerArray spacg_fn;

IntegerArray sixf_tp;
IntegerArray sixf_fp;

/********************/


void loadReadHeader(const std::string & fname,StringMap & Query, Integer2D & read_query,Integer2D & pos_read_query, IntegerMap & readlen);
void loadBWA(const std::string & fname, StringMap & Query, Integer2D & read_query,Integer2D & pos_read_query, IntegerMap & readlen, IntegerArray & query_mapped);
void loadBWAContig(const std::string & fname, StringMap & Query, Integer2D & read_query,Integer2D & pos_read_query, IntegerMap & readlen);
void loadFGSRef(std::string & fname, Integer2D & reads_fgs, IntegerArray & reads_all_fgs, StringMap & Query, const Integer2D & read_query, const Integer2D & pos_read_query,IntegerMap & readlen);
void loadFGS(std::string & fname, Integer2D & reads_fgs, Integer2D & reads_all_fgs, IntegerArray & reads_all_array, StringMap & Query,const Integer2D & read_query, const Integer2D & pos_read_query,IntegerMap & readlen, MapSet & pred_reads_map);
void loadFGSRead(std::string & fname, Integer2D & reads_fgs, Integer2D & reads_all_fgs, IntegerArray & reads_all_array, IntegerMap & readlen);
void evaluation(std::string & rocname, std::string & fname,  IntegerArray & actual_reads_fgs, Integer2D & all_pred_reads, Integer2D & pred_reads_fgs, IntegerArray & pred_reads_array, IntegerArray &all_tp, IntegerArray& all_fp, IntegerArray &all_fn);
void evaluation( IntegerArray & actual_reads_fgs, IntegerArray & all_pred_reads, IntegerArray & path_tp, IntegerArray & path_fp, IntegerArray & sixf_tp, IntegerArray & sixf_fp);
void getUniqueCount(IntegerArray &read_tp, IntegerArray &read_fp, IntegerArray &read_fn, IntegerArray &edge_tp, IntegerArray &edge_fp,IntegerArray &edge_fn , IntegerArray &path_tp, IntegerArray &path_fp,IntegerArray &path_fn, MapSet & pred_reads_map, IntegerArray & sixf_tp, IntegerArray & sixf_fp );
void addPerReadEdge(Integer2D & all_reads_path, Integer2D & all_reads, Integer2D & all_reads_edge);
void addPerReadEdge(Integer2D & reads_paths, Integer2D & reads_fgs, Integer2D & reads_edge);
void load_6f(const std::string & fname, IntegerArray & sixf_reads, Integer2D & all_reads, Integer2D & reads_6f, IntegerArray & pred_reads_array);
int ReadtoNumber(std::string & str, StringIntegerMap & NameMap);
int Name2Digit(std::string & str, StringMap & RefMap);
void sortPos(StringMap Query, Integer2D & read_query,Integer2D & pos_read_query, IntegerArray & query_mapped);
int upperIndex(const IntegerArray &arr, int n, int &y);
int lowerIndex(const IntegerArray &arr, int n, int &x);
