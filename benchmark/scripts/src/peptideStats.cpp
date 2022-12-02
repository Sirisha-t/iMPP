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
#include<boost/bimap.hpp>

typedef boost::bimap<std::string, int> StringIntMap;
typedef std::vector<std::string> StringArray;
typedef std::unordered_map<std::string, int> StringIntegerMap;
typedef std::map<int, int> IntIntMap;
typedef std::set<std::string> StringSet;
typedef std::set<int> IntegerSet;
typedef std::map<std::string, IntegerSet> String2DMap;
typedef std::vector< std::vector<int> > Integer2D;
typedef std::vector< std::vector<std::string> > String2D;
typedef std::vector<int> IntegerArray;


    StringIntegerMap alignedContLen60;
    StringIntegerMap alignedContLen70;
    StringIntegerMap alignedContLen80;
    StringIntegerMap alignedContLen90;

    StringIntegerMap alignedRefLen60;
    StringIntegerMap alignedRefLen70;
    StringIntegerMap alignedRefLen80;
    StringIntegerMap alignedRefLen90;


int NametoNumber(std::string & str, StringIntegerMap & NameMap)
{
        if(NameMap.count(str)==0)
                NameMap[str]=NameMap.size();

    return NameMap[str];
}

void RefLenMap(std::string & str, int len, StringIntegerMap & RefMap)
{
        if(RefMap.count(str)==0)
                RefMap[str]=len;
        else if(RefMap[str] < len)
                RefMap[str]= len;
}


void load_dmd(const std::string & refread,\
              const std::string & dmdin,\
              const std::string & contin,\
              StringIntegerMap & NameMap,\
              const std::string & dmdprot,\
              const std::string & outf,\
              const std::string & refin,\
              const std::string & contout,\
              int & cov)
{

    //IntegerSet Uniqread;
    StringSet Uniqread_60;
    StringSet Uniqread_70;
    StringSet Uniqread_80;
    StringSet Uniqread_90;
    StringIntegerMap contLength;
    StringIntegerMap refLength;
    String2DMap ReadsperContig;
    IntegerSet allReads;

    String2DMap ReadsperRef;
    IntegerSet allReadsinRef;




    IntegerSet readCount60;
    IntegerSet readCount70;
    IntegerSet readCount80;
    IntegerSet readCount90;
    IntegerSet refCount60 ;
    IntegerSet refCount70 ;
    IntegerSet refCount80 ;
    IntegerSet refCount90 ;
    int align_cont60 = 0;
    int align_cont70 = 0;
    int align_cont80= 0;
    int align_cont90 = 0;

    int align_ref60 = 0;
    int align_ref70 = 0;
    int align_ref80= 0;
    int align_ref90 = 0;

    StringSet ref60;
    StringSet ref70;
    StringSet ref80;
    StringSet ref90;
    //int ref60 = 0;
    //int ref70 = 0;
    //int ref80 = 0;
    //int ref90 = 0;
    int contig_count=0;
    int readcount = 0;
    int refcount = 0;
    int allreads_count = 0;
    int total_len = 0;
    long int total_ref_len = 0;
    //StringSet Uniqread;
    std::ifstream dxin;
    int discard=0;

        std::ifstream drin;
        drin.open(refread.c_str());
        if(!drin.is_open())
            std::cout<<"Cannot open file\n";
        else{
          //int count=0;

          int rnum=0;
          std::string line;
          boost::char_separator<char> sep("\t");

            while(getline(drin,line))
            {
               if(line[0]!='#')
               {

                 boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
                 boost::tokenizer<boost::char_separator<char> >::iterator it;
                 it = tokens.begin();
                 std::string read = (*it);
                  ++it;
                 std::string ref = *it;
                  ++it;
                   double seqid = std::stod(*it);
                    if(seqid >= 90.00)
                       {
                          /* if(NameMap.count(contig) == 0)
                              NameMap[contig] = 1;
                           else
                           {
                               int count = NameMap[contig];
                                  ++count;
                                     NameMap[contig] = count;
                           }*/
                           int rnum = NametoNumber(read, NameMap);
                           //ReadsperRef[ref].insert(rnum);
                           allReadsinRef.insert(rnum);

                       }
              }
          }
      }
	std::cout<<"Completed reading ref-read dmd\n";
       allreads_count = allReadsinRef.size();
    dxin.open(dmdin.c_str());
    if(!dxin.is_open())
        std::cout<<"Cannot open file\n";
    else{
      //int count=0;
      int discard=0;
      int rnum=0;
      std::string line;
      boost::char_separator<char> sep("\t");

        while(getline(dxin,line))
        {
           if(line[0]!='#')
           {

             boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
             boost::tokenizer<boost::char_separator<char> >::iterator it;
             it = tokens.begin();
             std::string read = (*it);
	            ++it;
             std::string contig = *it;
	            ++it;
	             double seqid = std::stod(*it);
	              if(seqid >= 90.00)
	                 {
		                  /* if(NameMap.count(contig) == 0)
			                    NameMap[contig] = 1;
		                   else
		                   {
		                       int count = NameMap[contig];
			                        ++count;
			                           NameMap[contig] = count;
		                   }*/
                       int rnum = NametoNumber(read, NameMap);
                       ReadsperContig[contig].insert(rnum);

	                 }
	        }
	    }
	}

    std::cout<<" Completed reading read as query file \n";
    std::ifstream fin, fin2, fin3;
    fin2.open(contin.c_str());
    std::ofstream fout, fout2;
    fout.open(outf.c_str());
    fout2.open(contout.c_str());
    if(!fin2.is_open())
        std::cout<<"Cannot open file\n";
    else{

      std::string line;

      //boost::char_separator<char> sep("\t");

        while(getline(fin2,line))
        {
           if(line[0]=='>')
           {
	            line.erase(0,1);
              std::string cont = line;
	            getline(fin2,line);
	            int len = line.length();
	            if(len>=cov)
              {
                //totalReads.insert(ReadsperContig[cont]);

                  for(auto it =  ReadsperContig[cont].begin(); it!= ReadsperContig[cont].end(); ++it)
                  {
                    allReads.insert(*it);
                  }

               total_len += len;
               contig_count++;
               //fout2<<">"<<cont<<"\n";
               //fout2<<line<<"\n";
             }
	            if(contLength.count(cont) == 0)
		              contLength[cont] = len;

	          }
        }
         readcount = allReads.size();
	       fout<<"# contigs (>="<<cov<<"nt) : "<<contig_count<<"\n";
         fout<<"# assembled reads(>="<<cov<<"nt) : "<<readcount<<"\n";
         fout<<"# total length : "<<total_len<<"\n";
     }
     std::cout<<"Completed storing contig lengths\n";

     fin3.open(refin.c_str());
     if(!fin3.is_open())
         std::cout<<"Cannot open file\n";
     else{
       int count=0;
       std::string line;

       //boost::char_separator<char> sep("\t");

         while(getline(fin3,line))
         {
            if(line[0]=='>')
            {
              refcount++;
 	            line.erase(0,1);
              std::string ref = line;
 	            getline(fin3,line);
 	            int len = line.length();
 	            if(refLength.count(ref) == 0)
 		              refLength[ref] = len;
              total_ref_len += len;

 	          }
         }
      }
      std::cout<<"Completed storing ref lengths\n";
      std::cout<<"Total predicted ref length : "<<total_ref_len<<"\n";

    fin.open(dmdprot.c_str());
    if(!fin.is_open())
        std::cout<<"Cannot open file\n";
    else{
      int count=0;
      int discard=0;
      int rnum=0;
      std::string line;
      boost::char_separator<char> sep("\t");

       	while(getline(fin,line))
        {
           if(line[0]!='#')
           {
             count = NameMap.size()+1;
             boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
             boost::tokenizer<boost::char_separator<char> >::iterator it;
             it = tokens.begin();
	            std::string query = (*it);
     	        ++it;
              std::string ref = (*it);
	            ++it;
             double percent = std::stod(*it);
	            ++it;
	            ++it;
	            ++it;
              ++it;
	            int q_start = std::stoi(*it);
	            ++it;
	            int q_end = std::stoi(*it);
              ++it;
              int r_start = std::stoi(*it);
              ++it;
              int r_end = std::stoi(*it);
	            int align_len = (abs(q_end - q_start + 1));
	            int actual_length = contLength[query];
              int ref_align_len = (abs(r_end - r_start + 1));
              int actual_ref_len = refLength[ref];
	            //double thresh = (0.6)*actual_length;

	          if( actual_length >= cov)
	           {
              if(align_len >= (0.40)*actual_length  )
                    {Uniqread_60.insert(query);
                    //ref70.insert(ref);
                    RefLenMap(ref, ref_align_len, alignedRefLen60);
                    RefLenMap(query, align_len, alignedContLen60);
                    }
	             if(align_len >= (0.50)*actual_length)
	              		{Uniqread_70.insert(query);
                      RefLenMap(ref, ref_align_len, alignedRefLen70);
                      RefLenMap(query, align_len, alignedContLen70);
                    //ref70.insert(ref);
                  }
	             if(align_len >= (0.60)*actual_length)
		             {Uniqread_80.insert(query);
                RefLenMap(ref, ref_align_len, alignedRefLen80);
                RefLenMap(query, align_len, alignedContLen80);
                 //ref70.insert(ref);

             }
	             if(align_len >= (0.70)*actual_length)
		             {Uniqread_90.insert(query);
                   RefLenMap(ref, ref_align_len, alignedRefLen90);
                   RefLenMap(query, align_len, alignedContLen90);
                 //ref70.insert(ref);
               }

              }
              // for ref seq match
            /*  if(ref_align_len >= (0.6)*actual_ref_len)
                  ref60 += 1;
              if(ref_align_len >= (0.7)*actual_ref_len)
                  ref70 += 1;
              if(ref_align_len >= (0.8)*actual_ref_len)
                  ref80 += 1;
              if(ref_align_len >= (0.9)*actual_ref_len)
                  ref90 += 1;*/

        }
      }

      for(StringIntegerMap::iterator i=alignedContLen60.begin(); i != alignedContLen60.end(); i++)
      {
        align_cont60 += i->second;
      }
      for(StringIntegerMap::iterator i=alignedContLen70.begin(); i != alignedContLen70.end(); i++)
      {
        align_cont70 += i->second;
      }
      for(StringIntegerMap::iterator i=alignedContLen80.begin(); i != alignedContLen80.end(); i++)
      {
        align_cont80 += i->second;
      }
      for(StringIntegerMap::iterator i=alignedContLen90.begin(); i != alignedContLen90.end(); i++)
      {
        align_cont90 += i->second;
      }


            for(StringIntegerMap::iterator i=alignedRefLen60.begin(); i != alignedRefLen60.end(); i++)
            {
              align_ref60 += i->second;
            }
            for(StringIntegerMap::iterator i=alignedRefLen70.begin(); i != alignedRefLen70.end(); i++)
            {
              align_ref70 += i->second;
            }
            for(StringIntegerMap::iterator i=alignedRefLen80.begin(); i != alignedRefLen80.end(); i++)
            {
              align_ref80 += i->second;
            }
            for(StringIntegerMap::iterator i=alignedRefLen90.begin(); i != alignedRefLen90.end(); i++)
            {
              align_ref90 += i->second;
            }

	for(StringSet::iterator it= Uniqread_60.begin(); it!= Uniqread_60.end();it++)
		{
      for(IntegerSet::iterator it2 = ReadsperContig[*it].begin(); it2 != ReadsperContig[*it].end(); ++it2)
        readCount60.insert(*it2);
      //readCount60 = readCount60 + NameMap[*it];
    }
	for(StringSet::iterator it= Uniqread_70.begin(); it!= Uniqread_70.end();it++)
  {
    for(IntegerSet::iterator it2 = ReadsperContig[*it].begin(); it2 != ReadsperContig[*it].end(); ++it2)
      readCount70.insert(*it2);
    //readCount60 = readCount60 + NameMap[*it];
  }

	for(StringSet::iterator it= Uniqread_80.begin(); it!= Uniqread_80.end();it++)
  {
    for(IntegerSet::iterator it2 = ReadsperContig[*it].begin(); it2 != ReadsperContig[*it].end(); ++it2)
      readCount80.insert(*it2);
    //readCount60 = readCount60 + NameMap[*it];
  }
	for(StringSet::iterator it= Uniqread_90.begin(); it!= Uniqread_90.end();it++)
  {
    for(IntegerSet::iterator it2 = ReadsperContig[*it].begin(); it2 != ReadsperContig[*it].end(); ++it2)
      readCount90.insert(*it2);
    //readCount60 = readCount60 + NameMap[*it];
  }

  /* for(StringSet::iterator it= ref60.begin(); it!= ref60.end();it++)
    {
      for(IntegerSet::iterator it2 = ReadsperRef[*it].begin(); it2 != ReadsperRef[*it].end(); ++it2)
        refCount60.insert(*it2);
      //readCount60 = readCount60 + NameMap[*it];
    }
    for(StringSet::iterator it= ref70.begin(); it!= ref70.end();it++)
      {
        for(IntegerSet::iterator it2 = ReadsperRef[*it].begin(); it2 != ReadsperRef[*it].end(); ++it2)
          refCount70.insert(*it2);
        //readCount60 = readCount60 + NameMap[*it];
      }

      for(StringSet::iterator it= ref80.begin(); it!= ref80.end();it++)
        {
          for(IntegerSet::iterator it2 = ReadsperRef[*it].begin(); it2 != ReadsperRef[*it].end(); ++it2)
            refCount80.insert(*it2);
          //readCount60 = readCount60 + NameMap[*it];
        }
        for(StringSet::iterator it= ref90.begin(); it!= ref90.end();it++)
          {
            for(IntegerSet::iterator it2 = ReadsperRef[*it].begin(); it2 != ReadsperRef[*it].end(); ++it2)
              refCount90.insert(*it2);
            //readCount60 = readCount60 + NameMap[*it];
          }
*/

      fin.clear();
      fin.close();
      fout<<"Seqcov"<<"\t"<<"Contig"<<"\t"<<"Aligned Ref"<<"\t"<<"Ref"<<"\t"<<"Cont-Tot"<<"\t"<<"Read"<<"\t"<<"All Read(Cont)"<<"\t"<<"All Read (Ref)"<<"\t"<<"Sp(cont)"<<"\t"<<"Sn(cont)"<<"\t"<<"Sp(read)"<<"\t"<<"Sn(read)\n";
      fout<<60<<"\t"<<align_cont60<<"\t"<<align_ref60<<"\t"<<total_ref_len<<"\t"<<total_len<<"\t"<<readCount60.size()<<"\t"<<readcount<<"\t"<<allreads_count<<"\t"<<((double)align_cont60/total_len)*100<<"\t"<<((double)align_ref60/total_ref_len)*100<<"\t"<<((double)readCount60.size()/readcount)*100<<"\t"<<((double)readCount60.size()/allreads_count)*100<<"\n" ;
      fout<<70<<"\t"<<align_cont70<<"\t"<<align_ref70<<"\t"<<total_ref_len<<"\t"<<total_len<<"\t"<<readCount70.size()<<"\t"<<readcount<<"\t"<<allreads_count<<"\t"<<((double)align_cont70/total_len)*100<<"\t"<<((double)align_ref70/total_ref_len)*100<<"\t"<<((double)readCount70.size()/readcount)*100<<"\t"<<((double)readCount70.size()/allreads_count)*100<<"\n";
	    fout<<80<<"\t"<<align_cont80<<"\t"<<align_ref80<<"\t"<<total_ref_len<<"\t"<<total_len<<"\t"<<readCount80.size()<<"\t"<<readcount<<"\t"<<allreads_count<<"\t"<<((double)align_cont80/total_len)*100<<"\t"<<((double)align_ref80/total_ref_len)*100<<"\t"<<((double)readCount80.size()/readcount)*100<<"\t"<<((double)readCount80.size()/allreads_count)*100<<"\n";
	    fout<<90<<"\t"<<align_cont90<<"\t"<<align_ref90<<"\t"<<total_ref_len<<"\t"<<total_len<<"\t"<<readCount90.size()<<"\t"<<readcount<<"\t"<<allreads_count<<"\t"<<((double)align_cont90/total_len)*100<<"\t"<<((double)align_ref90/total_ref_len)*100<<"\t"<<((double)readCount90.size()/readcount)*100<<"\t"<<((double)readCount90.size()/allreads_count)*100<<"\n";

      std::cout<<"File : "<<dmdin<<"\n";
      std::cout<<"Seqcov"<<"\t"<<"Contig"<<"\t"<<"Aligned Ref"<<"\t"<<"Ref"<<"\t"<<"Cont-Tot"<<"\t"<<"Read"<<"\t"<<"All Read(Cont)"<<"\t"<<"All Read (Ref)"<<"\t"<<"Sp(cont)"<<"\t"<<"Sn(cont)"<<"\t"<<"Sp(read)"<<"\t"<<"Sn(read)\n";
      std::cout<<60<<"\t"<<align_cont60<<"\t"<<align_ref60<<"\t"<<total_ref_len<<"\t"<<total_len<<"\t"<<readCount60.size()<<"\t"<<readcount<<"\t"<<allreads_count<<"\t"<<((double)align_cont60/total_len)*100<<"\t"<<((double)align_ref60/total_ref_len)*100<<"\t"<<((double)readCount60.size()/readcount)*100<<"\t"<<((double)readCount60.size()/allreads_count)*100<<"\n" ;
      std::cout<<70<<"\t"<<align_cont70<<"\t"<<align_ref70<<"\t"<<total_ref_len<<"\t"<<total_len<<"\t"<<readCount70.size()<<"\t"<<readcount<<"\t"<<allreads_count<<"\t"<<((double)align_cont70/total_len)*100<<"\t"<<((double)align_ref70/total_ref_len)*100<<"\t"<<((double)readCount70.size()/readcount)*100<<"\t"<<((double)readCount70.size()/allreads_count)*100<<"\n";
	    std::cout<<80<<"\t"<<align_cont80<<"\t"<<align_ref80<<"\t"<<total_ref_len<<"\t"<<total_len<<"\t"<<readCount80.size()<<"\t"<<readcount<<"\t"<<allreads_count<<"\t"<<((double)align_cont80/total_len)*100<<"\t"<<((double)align_ref80/total_ref_len)*100<<"\t"<<((double)readCount80.size()/readcount)*100<<"\t"<<((double)readCount80.size()/allreads_count)*100<<"\n";
	    std::cout<<90<<"\t"<<align_cont90<<"\t"<<align_ref90<<"\t"<<total_ref_len<<"\t"<<total_len<<"\t"<<readCount90.size()<<"\t"<<readcount<<"\t"<<allreads_count<<"\t"<<((double)align_cont90/total_len)*100<<"\t"<<((double)align_ref90/total_ref_len)*100<<"\t"<<((double)readCount90.size()/readcount)*100<<"\t"<<((double)readCount90.size()/allreads_count)*100<<"\n";

}
}

int main(int argc, char * argv[])
{
 	std::string refin = argv[1];
	std::string dmd_in = argv[2]; //output of diamond blastx (reads as query)
	std::string dmd_protdb = argv[3]; //diamond mapping (blastp) of contigs against uniprot db
	std::string contig_len_in = argv[4]; // faa file (for length of contig)
  	std::string ref_prot_len = argv[5]; // faa file for ref protein length
  	std::string coverage = argv[6];
  	int cov = std::stoi(coverage);
	std::string outfile = dmd_in.substr(0,dmd_in.find_last_of(".")) + ".seqlen." + coverage + ".count";
  	std::string contig_30 = contig_len_in.substr(0,contig_len_in.find_last_of(".")) + "." + coverage + ".fas";

	//StringIntMap ReadNameMap;
	StringIntegerMap ReadNameMap;
 	 //StringIntegerMap RefMap;
	std::cout<<"Reading file..\n";
	load_dmd(refin, dmd_in,contig_len_in,ReadNameMap,dmd_protdb, outfile, ref_prot_len, contig_30, cov);

	return 0;
}
