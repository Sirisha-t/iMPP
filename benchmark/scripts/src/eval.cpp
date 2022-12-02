#include "eval.h"

void evaluation( IntegerArray & actual_reads_fgs, IntegerArray & all_pred_reads, IntegerArray & path_tp, IntegerArray & path_fp, IntegerArray & sixf_tp, IntegerArray & sixf_fp)
{

  int tp=0;
  int fp=0;
  int fn=0;
  int tp_6 = 0;
  int fp_6=0;
  //std::cout<<"tp\tfp\tfn\tprecision\trecall\t1-prec\taccuracy\n";
  //std::cout<<"0\t0\t0\t0\t0\t0\t0\n";
  std::sort(actual_reads_fgs.begin(), actual_reads_fgs.end());
  std::sort(all_pred_reads.begin(), all_pred_reads.end());
  std::sort(path_tp.begin(), path_tp.end());
  std::sort(path_fp.begin(), path_fp.end());

  actual_reads_fgs.erase(std::unique(actual_reads_fgs.begin(),actual_reads_fgs.end()),actual_reads_fgs.end() );
  all_pred_reads.erase(std::unique(all_pred_reads.begin(),all_pred_reads.end()),all_pred_reads.end() );

  path_tp.erase(std::unique(path_tp.begin(),path_tp.end()),path_tp.end() );
  path_fp.erase(std::unique(path_fp.begin(),path_fp.end()),path_fp.end() );

  IntegerArray temp;
  IntegerArray temp2;
  IntegerArray all_tp;

  /***************************True Positives***************************************/
  std::set_intersection(path_tp.begin(),\
                          path_tp.end(),\
                          all_pred_reads.begin(),\
                          all_pred_reads.end(),\
                          std::back_inserter(sixf_tp));

  //std::set_difference( path_tp.begin(), path_tp.end(), temp.begin(), temp.end(),std::back_inserter(all_tp));
  //tp = all_tp.size();
  tp_6 = sixf_tp.size();
  //for(int i=0;i<all_tp.size(); i++)
    //std::cout<<all_tp[i]<<"\t"<<readNames.right.at(all_tp[i])<<"\n";
  temp.clear();

  /**************************False Positives****************************************/
  std::set_intersection(all_pred_reads.begin(),\
                    all_pred_reads.end(),\
                    path_fp.begin(),\
                    path_fp.end(),\
                    std::back_inserter(sixf_fp));
  //std::set_difference( path_fp.begin(), path_fp.end() ,temp.begin(), temp.end(),std::back_inserter(temp2));
  //fp = temp2.size();
  fp_6 = sixf_fp.size();
  temp.clear();
  temp2.clear();
  /**************************False Negatives*****************************************/
  std::set_difference(actual_reads_fgs.begin(),\
                    actual_reads_fgs.end(),\
                    all_pred_reads.begin(),\
                    all_pred_reads.end(),\
                    std::back_inserter(temp));
  std::set_difference(actual_reads_fgs.begin(), actual_reads_fgs.end(), all_tp.begin(), all_tp.end(), std::back_inserter(temp2));
  fn = temp2.size();
  temp.clear();
  //std::cout<<"six frame:\t"<<tp_6<<"\t"<<fp_6<<"\n";
  //std::cout<<"only path:\t"<<tp<<"\t"<<fp<<"\n";
  //std::cout<<"Computed matrix successfully!\n";

}

//Method to convert store read as integer value
int ReadtoNumber(std::string & str, StringIntegerMap & NameMap)
{
	int count = NameMap.size()+1;
	StringIntegerMap::left_const_iterator it = NameMap.left.find(str);
  	if(it == NameMap.left.end())
  		NameMap.insert(StringIntegerMap::value_type(str, count));
		//NameMap.left[str] =  NameMap.size();

    return NameMap.left.at(str);
}

void getUniqueCount(IntegerArray &read_tp, IntegerArray &read_fp, IntegerArray &read_fn, IntegerArray &edge_tp, IntegerArray &edge_fp,IntegerArray &edge_fn , IntegerArray &path_tp, IntegerArray &path_fp,IntegerArray &path_fn, MapSet & pred_reads_map, IntegerArray & sixf_tp, IntegerArray & sixf_fp )
{
  IntegerArray temp;
  IntegerArray tp, fp, fn, tpr, fpr;
  IntegerSet fp_path;
  std::sort(read_tp.begin(), read_tp.end());
  std::sort(read_fp.begin(), read_fp.end());

  read_tp.erase(std::unique(read_tp.begin(),read_tp.end()),read_tp.end() );
  read_fp.erase(std::unique(read_fp.begin(),read_fp.end()),read_fp.end() );

  std::sort(edge_tp.begin(), edge_tp.end());
  std::sort(edge_fp.begin(), edge_fp.end());

  edge_tp.erase(std::unique(edge_tp.begin(),edge_tp.end()),edge_tp.end() );
  edge_fp.erase(std::unique(edge_fp.begin(),edge_fp.end()),edge_fp.end() );

  std::sort(path_tp.begin(), path_tp.end());
  std::sort(path_fp.begin(), path_fp.end());

  path_tp.erase(std::unique(path_tp.begin(),path_tp.end()),path_tp.end() );
  path_fp.erase(std::unique(path_fp.begin(),path_fp.end()),path_fp.end() );

  std::cout<<"\t"<<"tp\tfp\tfn\n";
  /* getting unique edge count*/
  std::set_difference(edge_tp.begin(),\
                    edge_tp.end(),\
                    read_tp.begin(),\
                    read_tp.end(),\
                    std::back_inserter(temp));
  std::cout<<"edge:\t"<<temp.size()<<"\t";
  temp.clear();
  std::set_difference(edge_fp.begin(),\
                    edge_fp.end(),\
                    read_fp.begin(),\
                    read_fp.end(),\
                    std::back_inserter(temp));
  std::cout<<temp.size()<<"\t";
  temp.clear();
  std::set_difference(edge_fn.begin(),\
                    edge_fn.end(),\
                    read_fn.begin(),\
                    read_fn.end(),\
                    std::back_inserter(temp));
  std::cout<<temp.size()<<"\n";
  temp.clear();
  /************************/
  /* getting unique path count*/
  std::set_union(edge_tp.begin(),\
                    edge_tp.end(),\
                    read_tp.begin(),\
                    read_tp.end(),\
                    std::back_inserter(tp));

  std::set_union(edge_fp.begin(),\
                    edge_fp.end(),\
                    read_fp.begin(),\
                    read_fp.end(),\
                    std::back_inserter(fp));

  std::set_union(edge_fn.begin(),\
                    edge_fn.end(),\
                    read_fn.begin(),\
                    read_fn.end(),\
                    std::back_inserter(fn));


  std::set_union(sixf_tp.begin(),\
                 sixf_tp.end(),\
                 tp.begin(),\
                 tp.end(),\
                 std::back_inserter(tpr));

  std::set_union(sixf_fp.begin(),\
                sixf_fp.end(),\
                fp.begin(),\
                fp.end(),\
                std::back_inserter(fpr));

  std::sort(tpr.begin(), tpr.end());
  std::sort(fpr.begin(), fpr.end());
  tpr.erase(std::unique(tpr.begin(),tpr.end()),tpr.end() );
  fpr.erase(std::unique(fpr.begin(),fpr.end()),fpr.end() );
  /****/
  std::set_difference(path_tp.begin(),\
                      path_tp.end(),\
                      tpr.begin(),\
                      tpr.end(),\
                      std::back_inserter(temp));
  std::cout<<"path:\t"<<temp.size()<<"\t";
  temp.clear();
  std::set_difference(path_fp.begin(),\
                      path_fp.end(),\
                      fpr.begin(),\
                      fpr.end(),\
                      std::back_inserter(temp));
  std::cout<<temp.size()<<"\t";
  /*for(int i=0;i<temp.size(); i++)
  {
    for(const auto &it : pred_reads_map )
    {
      if(it.first == temp[i])
      {
        std::cout<<it.first<<" : "<<"\t";
        for(const auto & elem : it.second )
        {
            std::cout<<elem<<"\t";
            fp_path.insert(elem);
        }
        std::cout<<"\n";
      }
    }
  }*/
  //std::cout<<"fp paths size :"<<fp_path.size()<<"\n";
  //for(const auto &fp : fp_path)
    //std::cout<<fp<<"\t";
  std::cout<<"\n";

  std::cout<<"six fr:\t"<<sixf_tp.size()<<"\t"<<sixf_fp.size()<<"\n";
  temp.clear();
  /*std::set_difference(path_fn.begin(),\
                      path_fn.end(),\
                      fn.begin(),\
                      fn.end(),\
                      std::back_inserter(temp));
  std::cout<<temp.size()<<"\n";
  temp.clear();*/
  /************************/

}


int upperIndex(const IntegerArray &arr, int n, int &y)
{
    int l = 0, h = n - 1;
    while (l <= h) {
        int mid = (l + h) / 2;
        if (arr[mid] <= y)
            l = mid + 1;
        else
            h = mid - 1;
    }
    return h;
}
int lowerIndex(const IntegerArray &arr, int n, int &x)
{
    int l = 0, h = n - 1;
    while (l <= h) {
        int mid = (l + h) / 2;
        if (arr[mid] >= x)
            h = mid - 1;
        else
            l = mid + 1;
    }
    return l;
}


void evaluation(std::string & rocname, std::string & fname,  IntegerArray & actual_reads_fgs, Integer2D & all_pred_reads, Integer2D & pred_reads_fgs, IntegerArray & pred_reads_array, IntegerArray  &all_tp, IntegerArray& all_fp, IntegerArray &all_fn)
{
  std::ofstream fout1(fname.c_str());
  std::ofstream fout2(rocname.c_str());
  int tp=0;
  int fp=0;
  int fn=0;

  int tpr=0;
  int fpr=0;
  int fnr=0;
  //tp_out.resize(freq); // for ROC
  //tp_out.resize(1); // for benchmark
  fout1<<"tp\tfp\tfn\tprecision\trecall\t1-prec\taccuracy\n";
  fout1<<"0\t0\t0\t0\t0\t0\t0\n";
  fout2<<"tp\tfp\tfn\tprecision\trecall\t1-prec\taccuracy\n";
  fout2<<"0\t0\t0\t0\t0\t0\t0\n";
  std::sort(actual_reads_fgs.begin(), actual_reads_fgs.end());
  std::sort(pred_reads_array.begin(), pred_reads_array.end());

  actual_reads_fgs.erase(std::unique(actual_reads_fgs.begin(),actual_reads_fgs.end()),actual_reads_fgs.end() );
  pred_reads_array.erase(std::unique(pred_reads_array.begin(),pred_reads_array.end()),pred_reads_array.end() );

  IntegerArray temp;
  //IntegerArray all_tp;
  //IntegerArray all_fp;
  //IntegerArray all_fn;

  /***************************True Positives***************************************/
  std::set_intersection(actual_reads_fgs.begin(),\
                          actual_reads_fgs.end(),\
                          pred_reads_array.begin(),\
                          pred_reads_array.end(),\
                          std::back_inserter(all_tp));
  /**************************False Positives****************************************/
  std::set_difference(pred_reads_array.begin(),\
                    pred_reads_array.end(),\
                    actual_reads_fgs.begin(),\
                    actual_reads_fgs.end(),\
                    std::back_inserter(all_fp));

  /**************************False Negatives*****************************************/
  std::set_difference(actual_reads_fgs.begin(),\
                    actual_reads_fgs.end(),\
                    pred_reads_array.begin(),\
                    pred_reads_array.end(),\
                    std::back_inserter(all_fn));



  float prec = all_tp.size()/float(all_tp.size()+all_fp.size());
  float rec = all_tp.size()/float(all_tp.size()+all_fn.size());
  float inv_prec = 1 - float(prec);
  float acc = 2 * ( float(prec * rec)/ ( prec + rec));
  std::cout<<all_tp.size()<<"\t"<<all_fp.size()<<"\t"<<all_fn.size()<<"\t"<<prec<<"\t"<<rec<<"\t"<<acc<<"\n";
  for(int i=0;i<freq;i++){ //for ROC
  //for(int i=0;i<1;i++){ // for benchmark
  std::sort(pred_reads_fgs[i].begin(), pred_reads_fgs[i].end());
  std::sort(all_pred_reads[i].begin(), all_pred_reads[i].end());
  std::vector<int>::iterator it;
  std::vector<int>::iterator it2;
  pred_reads_fgs[i].erase(std::unique(pred_reads_fgs[i].begin(),pred_reads_fgs[i].end()),pred_reads_fgs[i].end() );
  all_pred_reads[i].erase(std::unique(all_pred_reads[i].begin(),all_pred_reads[i].end()),all_pred_reads[i].end() );

  std::set_intersection(actual_reads_fgs.begin(),\
                        actual_reads_fgs.end(),\
                        pred_reads_fgs[i].begin(),\
                        pred_reads_fgs[i].end(),\
                        std::back_inserter(temp));

  tp = temp.size();
  temp.clear();
  std::set_difference(
                      pred_reads_fgs[i].begin(),\
                      pred_reads_fgs[i].end(),\
                      actual_reads_fgs.begin(),\
                      actual_reads_fgs.end(),\
                      std::back_inserter(temp));

  fp = temp.size();
  temp.clear();
  std::set_difference(actual_reads_fgs.begin(),\
                  actual_reads_fgs.end(),\
                  pred_reads_fgs[i].begin(),\
                  pred_reads_fgs[i].end(),\
                  std::back_inserter(temp));
  fn = temp.size();
  temp.clear();

  std::set_intersection(actual_reads_fgs.begin(),\
                        actual_reads_fgs.end(),\
                        all_pred_reads[i].begin(),\
                        all_pred_reads[i].end(),\
                        std::back_inserter(temp));

  tpr = temp.size();
  temp.clear();
  std::set_difference(
                      all_pred_reads[i].begin(),\
                      all_pred_reads[i].end(),\
                      actual_reads_fgs.begin(),\
                      actual_reads_fgs.end(),\
                      std::back_inserter(temp));

  fpr = temp.size();
  temp.clear();
  std::set_difference(actual_reads_fgs.begin(),\
                  actual_reads_fgs.end(),\
                  all_pred_reads[i].begin(),\
                  all_pred_reads[i].end(),\
                  std::back_inserter(temp));
  fnr = temp.size();
  temp.clear();



   float prec = tp/float(tp+fp);
   float rec = tp/float(tp+fn);
   float inv_prec = 1 - float(prec);
   float acc = 2 * ( float(prec * rec)/ ( prec + rec));


   float precr = tpr/float(tpr+fp);
   float recr = tpr/float(tpr+fnr);
   float inv_precr = 1 - float(precr);
   float accr = 2 * ( float(precr * recr)/ ( precr + recr));

  fout1<<tp<<"\t"<<fp<<"\t"<<fn<<"\t"<<prec<<"\t"<<rec<<"\t"<<inv_prec<<"\t"<<acc<<"\n";
  fout2<<tpr<<"\t"<<fp<<"\t"<<fnr<<"\t"<<precr<<"\t"<<recr<<"\t"<<inv_precr<<"\t"<<accr<<"\n";

  }
  fout1<<"1\t1\t1\t1\t1\t1\t1\n";
  fout2<<"1\t1\t1\t1\t1\t1\t1\n";
  fout1.close();
  fout2.close();
  std::cout<<"Computed matrix successfully!\n";

}

int Name2Digit(std::string & str, StringMap & RefMap)
{
  	if(RefMap.count(str)==0)
  		RefMap[str]=RefMap.size();

    return RefMap[str];
}

void load_6f(const std::string & fname, IntegerArray & sixf_reads, Integer2D & all_reads, Integer2D & reads_6f, IntegerArray & pred_reads_array)
{

  std::ifstream fin;
  fin.open(fname.c_str());
  int check=0;
  if(!fin.is_open())
        {std::cout<<"Cannot open file!"<<"\n";
        exit(1);}
  else
  {
      std::string line;
      while(std::getline(fin,line))
      {
      int L=0;
      if(line[0] == '>')
      {
        line.erase(0,1);
        int read_name = readNames.left.at(line);
        sixf_reads.push_back(read_name);
        pred_reads_array.push_back(read_name);
        for(int k=0;k<freq;k++)
          {
            reads_6f[k].push_back(read_name);
            all_reads[k].push_back(read_name);
          }
      }
        else
        {
                std::string seq = line;

        }
    }
}
}


void addAllReadEdge(IntegerArray & all_reads_path, IntegerArray & all_reads, IntegerArray & all_reads_edge)
{

    for(int i=0; i<all_reads.size();i++)
    {
      all_reads_path.push_back(all_reads[i]);
    }
    for(int k=0;k<all_reads_edge.size();k++)
        all_reads_path.push_back(all_reads_edge[k]);

    for(int i=0; i<all_reads.size();i++)
    {
      all_reads_edge.push_back(all_reads[i]);
    }

}

void addPerReadEdge(Integer2D & reads_paths, Integer2D & reads_fgs, Integer2D & reads_edge)
{
  std::cout<<"Adding orphan reads to edges and paths..\n";
  //std::cout<<"Initial size of path "<<reads_paths.size()<<"\n";
  for(int i=0;i<reads_paths.size();i++)
  {
    for(int j=0; j<reads_fgs[i].size();j++)
    {
      reads_paths[i].push_back(reads_fgs[i][j]);
    }
    for(int k=0;k<reads_edge[i].size();k++)
        reads_paths[i].push_back(reads_edge[i][k]);
  }
  for(int i=0;i<reads_edge.size();i++)
  {
    for(int j=0; j<reads_fgs[i].size();j++)
    {
      reads_edge[i].push_back(reads_fgs[i][j]);
    }
  }


}

void loadFGS(std::string & fname,
  Integer2D & reads_fgs,
  Integer2D & reads_all_fgs,
  IntegerArray & reads_all_array,
  StringMap & Query,
  const Integer2D & read_query,
  const Integer2D & pos_read_query,
  IntegerMap & readlen,
  MapSet & pred_reads_map)
 {
  std::ifstream fin;
  int check=0;
  fin.open(fname.c_str());
  if(!fin.is_open())
      {
        std::cout<<"Cannot open file "<<fname<<"\n";
        exit(1);
      }
  else
  {
      boost::char_separator<char> sep("\t");
      std::cout<<"Reading fgs file "<<fname<<"\n\n";
      std::string line;
      check++;
      if(check==1)
        {
          reads_fgs.resize(freq); // for ROC
          reads_all_fgs.resize(freq);
        }
      if(std::getline(fin,line))
      {
        if(line[0]=='>')
        {

          line.erase(0,1);
          std::string query = line;
          int query_name=Name2Digit(query,Query);
          while(std::getline(fin,line))
          {
            //std::cout<<query<<"\t"<<query_name<<"\n";
            if( line[0] != '>'){
              //std::cout<<line<<"\n";
              boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
              boost::tokenizer<boost::char_separator<char> >::iterator it;
              it = tokens.begin();
              int start=std::stoi(*it);
              ++it;
              int end= std::stoi(*it);
              if(start>end)
              {
                  int temp=start;
                  start=end;
                  end=temp;
              }
              ++it;
              ++it;
              ++it;
              float score = std::stof(*it);

              for(int i=0;i<read_query[query_name].size();i++)
              {
                  int read_name=read_query[query_name][i];
                  int pos=pos_read_query[query_name][i];
                  int L = readlen[read_name];
                  int overlap = L*threshold + 1 ;
                  float minv = 1.2;
                  float minv_roc = 1.2;
                  float maxv = 1.4;
                  float maxv_roc = 1.5;
                  float minv_inc = 0.0;
                  float minv_inc_roc = 0.0;
                  float inc = (maxv - minv)/freq;
                  float inc_roc = (maxv_roc - minv_roc)/freq;
                  //double ch = double(end-start+1)/L;
                  int st = start-L+overlap;
		              int en = end-overlap+1;
		                /** for ROC with multiple values **/
                  //std::cout<<read_name<<"\t"<<pos<<"\t"<<L<<overlap<<"\n";
                    if((start-L+overlap <=pos) && (pos<= end-overlap+1))
                    {
                      //std::cout<<start<<"\t"<<end<<"\t"<<overlap<<"\t"<<pos<<"\n";
                      //if(score <= 1.395)

                        pred_reads_map[read_name].insert(query_name);
                        reads_all_array.push_back(read_name);
                        for(int i=0; i<freq; i++){
                         minv_inc = minv + inc;
                         minv_inc_roc = minv_roc + inc_roc;
                         if(score <= minv ){
                           reads_fgs[i].push_back(read_name);
                         }
                         if(score <= minv_roc ){
                           reads_all_fgs[i].push_back(read_name);
                         }
                         minv = minv + inc;
                         minv_roc = minv_roc + inc_roc;
                       }

                    }
              }
            }

            else{
              line.erase(0,1);
              query = line;
              query_name=Name2Digit(query,Query);
              }
            }
          }
        }
      }
      std::cout<<"Completed reading file "<<fname<<"\n\n";
    fin.clear();
    fin.close();
  }


void loadFGSRef(std::string & fname,
  Integer2D & reads_fgs,
  IntegerArray & reads_all_fgs,
  StringMap & Query,
  const Integer2D & read_query,
  const Integer2D & pos_read_query,
  IntegerMap & readlen )
 {
  std::ifstream fin;
  int check=0;
  int count =0;
  fin.open(fname.c_str());
  if(!fin.is_open()){
    std::cout<<"Cannot open file "<<fname<<"\n";
    exit(1);
  }

  else
  {
      boost::char_separator<char> sep("\t");
      std::cout<<"Reading fgs ref file "<<fname<<"\n\n";
      std::string line;
      check++;
      if(check==1)
        {
          reads_fgs.resize(freq); // for ROC

        }
      while(std::getline(fin,line))
      {
        if(line[0]!='#')
        {


              boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
              boost::tokenizer<boost::char_separator<char> >::iterator it;
              it = tokens.begin();
              std::string query = (*it);
              if(Query.count(query)!=0)
              {
                  int query_name=Name2Digit(query,Query);

              ++it;
              ++it;
              ++it;
              int start=std::stoi(*it);
              ++it;
              int end= std::stoi(*it);
              if(start>end)
              {
                int temp=start;
                start=end;
                end=temp;

              }

              int upper  = upperIndex(pos_read_query[query_name], pos_read_query[query_name].size(), end);
              int lower = lowerIndex(pos_read_query[query_name], pos_read_query[query_name].size(), start);
              for(int i = lower; i <= upper; i++)
              {
                 int read_name=read_query[query_name][i];
                  int pos=pos_read_query[query_name][i];
                  int L = readlen[read_name];
                  int overlap = L*threshold + 1 ;

                  if((start-L+overlap <=pos) && (pos<= end-overlap+1))
                    {
                      //if(score <= 1.395)
                        reads_all_fgs.push_back(read_name);
                        for(int j=0; j<freq; j++){
                           reads_fgs[j].push_back(read_name);
                         }

                    }
                  }

                }
              }
            }
          }

      std::cout<<"Completed reading "<<fname<<" file.\n\n";
    fin.clear();
    fin.close();
  }


void loadFGSRead(std::string & fname, Integer2D & reads_fgs, Integer2D & reads_all_fgs, IntegerArray & reads_all_array, IntegerMap & readlen)
{
    std::ifstream fin;
    fin.open(fname.c_str());
    int check=0;
    int L=0;
    if(!fin.is_open())
    	{
        std::cout<<"Cannot open "<<fname<<" file!"<<"\n";
        exit(1);
      }
    else
    {
        std::string line;
        boost::char_separator<char> sep("\t");
        std::cout<<"Reading fgs file "<<fname<<"\n\n";
        if(std::getline(fin,line))
        {
        if(line[0] == '>')
        {
          check++;
          if(check==1)
            {
              reads_fgs.resize(freq); // for ROC
              //reads_fgs.resize(1); // for benchmark
              reads_all_fgs.resize(freq);
            }
          line.erase(0,1);
          std::string query = line;
          int read_name=readNames.left.at(query);
          //std::cout<<query<<"\t"<<read_name<<"\n";
	         while(std::getline(fin,line))
          {
            if( line[0] != '>'){
              boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
              boost::tokenizer<boost::char_separator<char> >::iterator it;
              it = tokens.begin();
              int start=std::stoi(*it);
              ++it;
              int end=std::stoi(*it);
              if(readlen[read_name] > 0)
              {
      	       L=readlen[read_name];
                 if(start>end)
                 {
                   int temp=start;
                   start=end;
                   end=temp;
                 }
                 ++it;
                 ++it;
                 ++it;
                 float score = std::stof(*it);
                 float minv = 1.2;
                 float minv_roc = 1.2;
                 float maxv = 1.4;
                 float maxv_roc = 1.5;
                 float minv_inc = 0.0;
                 float minv_inc_roc = 0.0;
                 float inc = (maxv - minv)/freq;
                 float inc_roc = (maxv_roc - minv_roc)/freq;
                 double ch = double(end-start+1)/L;

                   /** for ROC **/
                   if(ch >= rthresh)
                   {
                     //if(score <= 1.395)
                     //std::cout<<read_name<<"\n";
                      reads_all_array.push_back(read_name);
                      for(int i=0; i<freq; i++){
                        minv_inc = minv + inc;
                        minv_inc_roc = minv_roc + inc_roc;
                        if(score <= minv ){
                          reads_fgs[i].push_back(read_name);
                        }
                        if(score <= minv_roc ){
                          reads_all_fgs[i].push_back(read_name);
                        }
                        minv = minv + inc;
                        minv_roc = minv_roc + inc_roc;
                      }
             }
            }

          }

            else{
              line.erase(0,1);
              query = line;
              read_name=readNames.left.at(query);
            }

      }
    }
  }
    std::cout<<"Completed reading "<<fname<<" file.\n\n";
    std::cout<<"Reads size : "<<reads_all_array.size()<<"\n";
    fin.clear();
    fin.close();
}
}


void sortPos(StringMap Query, Integer2D & read_query,Integer2D & pos_read_query)
{
   std::cout<<"Sorting pos\n";
    for(auto &x: Query)
    {
	std::cout<<"x:"<<x.second<<"\n";
        int query_size = read_query[x.second].size();
        std::cout<<"query size:"<<query_size<<"\n";
	std::pair<int, int> pairt[query_size];
        // Storing the respective array
        // elements in pairs.
        for (int i = 0; i < query_size; i++)
        {
        pairt[i].first = pos_read_query[x.second][i];
        pairt[i].second = read_query[x.second][i];
        }

        // Sorting the pair array.
        std::sort(pairt, pairt + query_size);

        // Modifying original arrays
        for (int i = 0; i < query_size; i++)
         {
            pos_read_query[x.second][i] = pairt[i].first;
            read_query[x.second][i] = pairt[i].second;
        }
    for(int i=0; i < query_size; i++)
    {
        std::cout<<pos_read_query[x.second][i]<<"\t"<<read_query[x.second][i]<<"\n";
    }


    }
}

void loadBWA(const std::string & fname,
  StringMap & Query,
  Integer2D & read_query,
  Integer2D & pos_read_query,
  IntegerMap & readlen)
{
    std::ifstream fin;
    fin.open(fname.c_str());
    if(!fin.is_open())
        {
          std::cout<<"Cannot open "<<fname<<" file\n";
          exit(1);
        }
    else{
      int count=0;
      int discard=0;
	    std::string line;
      boost::char_separator<char> sep("\t");
      std::cout<<"Reading bwa file "<<fname<<"\n\n";
      while(getline(fin,line)){
        if(line[0]=='@' && line[1]=='S'){
          boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
          boost::tokenizer<boost::char_separator<char> >::iterator it;
          it  = tokens.begin();
          ++it;
          std::string line = (*it);
          int colon = line.find(':');
          line.erase(0,colon+1);
          int query_num = Name2Digit(line, Query);
        }
        else if(line[0]!='@'){
          count++;
          if(count==1){
            std::cout<<"Query size:"<<Query.size()<<"\n";
            std::cout<<"Query names have been successfully stored\n\n";
            read_query.resize(Query.size()+1);
            pos_read_query.resize(Query.size()+1);
          }
          boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
          boost::tokenizer<boost::char_separator<char> >::iterator it;
          it = tokens.begin();
          std::string read = (*it);
          int read_name= ReadtoNumber(read,readNames);
	        ++it;
          int flag = std::stoi(*it);
         // if(flag != 256){
          ++it;
          if((*it)!= "*"){
            std::string query = (*it);
            int query_name=Name2Digit(query,Query);
            ++it;
            int position = std::stoi(*it);
            ++it;
            ++it;
            ++it;
            ++it;
            ++it;
            ++it;
            //readlen[read_name] = (*it).length();
            if(readlen.count(read_name)==0)
              readlen[read_name] = (*it).length();
            else if( readlen[read_name] < (*it).length())
                  readlen[read_name] = (*it).length();
            read_query[query_name].push_back(read_name);
            pos_read_query[query_name].push_back(position);
          //}
        }
      }
      else  discard++;
    }
  }
  std::cout<<"Completed reading "<<fname<<"file.\n\n";
  fin.clear();
  fin.close();
}

void loadBWAContig(const std::string & fname,
  StringMap & Query,
  Integer2D & read_query,
  Integer2D & pos_read_query,
  IntegerMap & readlen)
{
    std::ifstream fin;
    fin.open(fname.c_str());
    if(!fin.is_open())
        {
          std::cout<<"Cannot open "<<fname<<" file\n";
          exit(1);
        }
    else{
      int count=0;
      int discard=0;
	    std::string line;
      boost::char_separator<char> sep("\t");
      boost::char_separator<char> sep2(".");

      std::cout<<"Reading bwa file "<<fname<<"\n\n";
      while(getline(fin,line)){
        if(line[0]=='@' && line[1]=='S'){
          boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
          boost::tokenizer<boost::char_separator<char> >::iterator it;
          it  = tokens.begin();
          ++it;
          std::string line = (*it);
          int colon = line.find(':');
          line.erase(0,colon+1);
          int hyphen = line.find('_');
          if (line.substr(0,hyphen) == "NODE")
          {
            // spades contig file
            boost::tokenizer<boost::char_separator<char> > tokens(line, sep2);
            boost::tokenizer<boost::char_separator<char> >::iterator it;
            it  = tokens.begin();
            std::string node = *it;
            ++it;
            std::string cov = *it;
            int index = cov.find('_');
            cov.erase(index);
            std::string query = node + "."+cov;
            int query_num = Name2Digit(query, Query);
            //std::cout<<query<<"\t"<<line<<"\n";

          }
          else{
            // sga contig file
            line.erase(hyphen);
            //std::cout<<line<<"\n";
            int query_num = Name2Digit(line, Query);
          }

        }
        else if(line[0]!='@'){
          count++;
          int query_name;
          if(count==1){
            std::cout<<"Query size:"<<Query.size()<<"\n";
            std::cout<<"Query names have been successfully stored\n\n";
            read_query.resize(Query.size()+1);
            pos_read_query.resize(Query.size()+1);
          }
          boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
          boost::tokenizer<boost::char_separator<char> >::iterator it;
          it = tokens.begin();
          std::string read = (*it);
          int read_name= ReadtoNumber(read,readNames);
	        ++it;
          int flag = std::stoi(*it);
         // if(flag != 256){
          ++it;
          if((*it)!= "*"){
            std::string query = (*it);
            int hyphen = query.find('_');
            if (query.substr(0,hyphen) == "NODE")
              {
                // spades contig file
                boost::tokenizer<boost::char_separator<char> > tokens(query, sep2);
                boost::tokenizer<boost::char_separator<char> >::iterator it;
                it  = tokens.begin();
                std::string node = *it;
                ++it;
                std::string cov = *it;
                int index = cov.find('_');
                cov.erase(index);
                std::string query_sub = node + "."+cov;
                query_name = Name2Digit(query_sub, Query);
                //std::cout<<query_sub<<"\t"<<query_name<<"\n";

              }//std::cout<<line<<"\n";
            else{
              query.erase(hyphen);
              query_name = Name2Digit(query, Query);
            }

            //int query_name=Name2Digit(query,Query);
            ++it;
            int position = std::stoi(*it);
            ++it;
            ++it;
            ++it;
            ++it;
            ++it;
            ++it;
            //std::cout<<query<<"\t"<<query_name<<"\t"<<position<<"\t"<<read_name<<"\n";
                        //readlen[read_name] = (*it).length();
            if(readlen.count(read_name)==0)
              readlen[read_name] = (*it).length();
            else if( readlen[read_name] < (*it).length())
                  readlen[read_name] = (*it).length();
            read_query[query_name].push_back(read_name);
            pos_read_query[query_name].push_back(position);
          //}
        }
      }
      else  discard++;
    }
  }
  std::cout<<"Completed reading "<<fname<<"\n\n";
  fin.clear();
  fin.close();
}


bool readInputFileNames(std::ifstream & infile)
{
  std::string line;
  while(getline(infile, line))
  {
    line.erase(std::remove_if(line.begin(), line.end(), isspace),line.end());
    if (line[0] == '#' || line.empty()) continue;
    auto sep = line.find("=");
    auto name = line.substr(0, sep);
    auto value = line.substr(sep + 1);

    //storing file names
    if (name == "bwa_ref") bwa_ref = value;
    else if (name == "bwa_path") bwa_path = value;
    else if (name == "bwa_edge") bwa_edge = value;
    else if (name == "bwa_sgacg") bwa_sgacg = value;
    else if (name == "bwa_spadescg") bwa_spadescg = value;
    else if (name == "fgs_ref") fgs_ref = value;
    else if (name == "fgs_read") fgs_read = value;
    else if (name == "fgs_path") fgs_path = value;
    else if (name == "fgs_edge") fgs_edge = value;
    else if (name == "fgs_sgacg") fgs_sgacg = value;
    else if (name == "fgs_spadescg") fgs_spadescg = value;
    else if (name == "sixfr") sixfr_reads = value;
    else if (name == "fastq") reads_fastq = value;
}
  //Storing output eval file names
  eval_read = fgs_read.substr(0,fgs_read.find_last_of(".")) + ".eval.fgs.txt";
  eval_path = fgs_path.substr(0,fgs_path.find_last_of(".")) + ".eval.fgs.txt";
  eval_edge = fgs_edge.substr(0,fgs_edge.find_last_of(".")) + ".eval.fgs.txt";
  eval_sgacg = fgs_sgacg.substr(0,fgs_sgacg.find_last_of(".")) + ".eval.fgs.txt";
  eval_spadescg = fgs_spadescg.substr(0,fgs_spadescg.find_last_of(".")) + ".eval.fgs.txt";

  //storing roc output
  roc_read = fgs_read.substr(0,fgs_read.find_last_of(".")) + ".fgs.roc.txt";
  roc_path = fgs_path.substr(0,fgs_path.find_last_of(".")) + ".fgs.roc.txt";
  roc_edge = fgs_edge.substr(0,fgs_edge.find_last_of(".")) + ".fgs.roc.txt";
  roc_sgacg = fgs_sgacg.substr(0,fgs_sgacg.find_last_of(".")) + ".fgs.roc.txt";
  roc_spadescg = fgs_spadescg.substr(0,fgs_spadescg.find_last_of(".")) + ".fgs.roc.txt";

  infile.clear();
  infile.close();
  return true;
}


int main(int argc, char * argv[])
{
  if(argc < 2)
  {
      std::cout<<"ERROR:Input parameter file must be provided. \n";
  }
  else{
    std::ifstream infile;
    infile.open(argv[1]);
    if(!infile.is_open())
    {
      std::cout<<"ERROR:Cannot open input parameter file. Please check input.\n";
      exit(1);
    }
    else{
      if(readInputFileNames(infile)){
        /** reading ref files **/
        loadBWA(bwa_ref,refName,read_RefQuery,pos_of_readRef, lenRef);
        sortPos(refName,read_RefQuery,pos_of_readRef);
        loadFGSRef(fgs_ref,fgs_readsPerRef, fgs_allReadsRef, refName,read_RefQuery,pos_of_readRef, lenRef );


        /** eval for reads*/
        loadFGSRead(fgs_read, fgs_readsPerRead, fgs_allReads, fgs_allReadsArray, lenRef);
        evaluation(roc_read,eval_read,fgs_allReadsRef,fgs_allReads,fgs_readsPerRead, fgs_allReadsArray, read_tp, read_fp, read_fn);

        /** eval for sga **/
        loadBWAContig(bwa_sgacg,SGAcontigName,read_SGAContigQuery,pos_of_SGAreadContig, lenSGACG);
        sortPos(SGAcontigName,read_SGAContigQuery,pos_of_SGAreadContig);
        loadFGS(fgs_sgacg,fgs_readsPerSGACG, fgs_allReadsSGACG, fgs_allReadsSGAArray, SGAcontigName,read_SGAContigQuery,pos_of_SGAreadContig, lenSGACG, predReadsSGA);
        evaluation(roc_sgacg,eval_sgacg,fgs_allReadsRef,fgs_allReadsSGACG, fgs_readsPerSGACG, fgs_allReadsSGAArray, sgacg_tp, sgacg_fp, sgacg_fn );
        
        /** eval for spades **/
        loadBWAContig(bwa_spadescg,SPAcontigName,read_SPAContigQuery,pos_of_SPAreadContig, lenSPACG);
        sortPos(SPAcontigName,read_SPAContigQuery,pos_of_SPAreadContig);
        loadFGS(fgs_spadescg,fgs_readsPerSPACG, fgs_allReadsSPACG, fgs_allReadsSPAArray, SPAcontigName,read_SPAContigQuery,pos_of_SPAreadContig, lenSPACG, predReadsSPA);
        evaluation(roc_spadescg, eval_spadescg,fgs_allReadsRef,fgs_allReadsSPACG, fgs_readsPerSPACG, fgs_allReadsSPAArray, spacg_tp, spacg_fp, spacg_fn);


        /** eval for edges and paths */
        loadBWA(bwa_edge,edgeName,read_EdgeQuery,pos_of_readEdge, lenEdge);
        loadBWA(bwa_path,pathName,read_PathQuery,pos_of_readPath, lenPath);

        sortPos(edgeName,read_EdgeQuery,pos_of_readEdge);
        sortPos(pathName,read_PathQuery,pos_of_readPath);

        loadFGS(fgs_edge, fgs_readsPerEdge, fgs_allReadsEdge, fgs_allReadsEdgeArray, edgeName,read_EdgeQuery,pos_of_readEdge, lenEdge, predReadsEdge);
        loadFGS(fgs_path,fgs_readsPerPath, fgs_allReadsPath, fgs_allReadsPathArray, pathName,read_PathQuery,pos_of_readPath, lenPath, predReadsPath);
        addPerReadEdge(fgs_allReadsPath, fgs_allReads, fgs_allReadsEdge); // adding read and edge prediction to path(2D)
        addPerReadEdge(fgs_readsPerPath,fgs_readsPerRead,fgs_readsPerEdge); //adding  reads and edge prediction to path (1D)
        addAllReadEdge(fgs_allReadsPathArray, fgs_allReadsArray, fgs_allReadsEdgeArray);
        load_6f(sixfr_reads,sixf_reads,fgs_allReadsPath,fgs_readsPerPath, fgs_allReadsPathArray);
        evaluation(roc_edge,eval_edge,fgs_allReadsRef,fgs_allReadsEdge, fgs_readsPerEdge, fgs_allReadsEdgeArray, edge_tp, edge_fp, edge_fn);
        evaluation(roc_path,eval_path,fgs_allReadsRef,fgs_allReadsPath, fgs_readsPerPath, fgs_allReadsPathArray, path_tp, path_fp, path_fn);

        /** getting uniq count */
        evaluation(fgs_allReadsRef, sixf_reads, path_tp, path_fp, sixf_tp, sixf_fp);
        getUniqueCount(read_tp, read_fp, read_fn, edge_tp, edge_fp, edge_fn, path_tp, path_fp, path_fn, predReadsPath, sixf_tp, sixf_fp);



      }
    }
  }
  return 0;
}
