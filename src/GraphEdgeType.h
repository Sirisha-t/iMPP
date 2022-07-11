#ifndef __GRAPHEDGETYPE_H_
#define __GRAPHEDGETYPE_H_

#include <iostream>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <assert.h>

#include "DataType.h"   

class GraphEdgeType
{
  public:

    // empty construction function; setting the overlap length to 0
    GraphEdgeType() {
        overlap_ = cov_ = num_mismatch_ = 0;
        ov_pos_[0] = ov_pos_[1] = ov_pos_[2] = ov_pos_[3] = 0; 
        cigar_ = "";
        visited_ = resolved_ = is_rc_ = false;
        src_orientation_ = true;
    }

    // initialization with overlap length
    explicit GraphEdgeType(const SeqIdxType overlap)    {
        overlap_ = overlap;
        ov_pos_[0] = ov_pos_[1] = ov_pos_[2] = ov_pos_[3] = 0; 
        cov_ = num_mismatch_ = 0;
        cigar_ = "";
        visited_ = resolved_ = is_rc_ = false;
        src_orientation_ = true;
    }

    // destructor function
    ~GraphEdgeType()    {
        ;   // do nothing
    }

    // assignment operator
    GraphEdgeType& operator=(const GraphEdgeType &n) {
        this->overlap_ = n.overlap_;
        this->num_mismatch_ = n.num_mismatch_;
        memcpy(this->ov_pos_, n.ov_pos_, 4 * sizeof(SeqIdxType));
        this->cov_ = n.cov_;
        this->cigar_ = n.cigar_;
        this->src_orientation_ = n.src_orientation_;
        this->resolved_ = n.resolved_;
        this->visited_ = n.visited_;
        this->is_rc_ = n.is_rc_;
        return *this;
    } 

    // setting the overlap length
    void SetOverlap(const SeqIdxType overlap)   {
        overlap_ = overlap;
    }

    // retrieving the overlap length
    SeqIdxType GetOverlap(void) const {
        return overlap_;
    }

    // setting the coverage
    void SetCoverage(const CoverageType c)  {
        cov_ = c;
        return;
    }

    // returns the coverage of the current node
    CoverageType GetCoverage(void) const    {
        return cov_;
    }

    // setting the source orientation
    void SetSrcOrientation(const bool c)  {
        src_orientation_ = c;
        return;
    }

    // returns the coverage of the current node
    bool GetSrcOrientation(void) const    {
        return src_orientation_;
    }

    // returns whether the node has been visited or not
    bool IsVisited(void) const  {
        return visited_;
    }

    // setting the node's visited status
    void SetVisited(const bool b) {
        visited_ = b;
        return;
    }

    // returns whether the node has been resolved
    bool IsResolved(void) const  {
        return resolved_;
    }

    // setting the node's resolved status
    void SetResolved(const bool b) {
        resolved_ = b;
        return;
    }

    // setting the "is reverse complement" tag
    void SetIsRevComplement(const bool b)   {
        is_rc_ = b;
        return;
    }

    // returns whether the adjacent nodes overlap with reverse complementarity
    bool IsRevComplement(void) const    {
        return is_rc_;
    }

    // setting the num. of mismatches
    void SetNumMismatches(const int m)   {
        num_mismatch_ = m;
        return;
    }

    // returns the number of mismatches
    int GetNumMismatches(void) const    {
        return num_mismatch_;
    }

    // setting alignment CIGAR string
    void SetAlignmentCIGAR(const std::string & s)   {
        cigar_ = s;
        return;
    }

    // returns the alignment cigar string
    std::string GetAlignmentCIGAR(void) const    {
        return cigar_;
    }

    // setting overlapped positions
    void SetOverlapPosition(const SeqIdxType a, const SeqIdxType b, const SeqIdxType c, const SeqIdxType d)   {
        ov_pos_[0] = a;
        ov_pos_[1] = b;
        ov_pos_[2] = c;
        ov_pos_[3] = d;
        return;
    }

    // returns the ith value of the overlap position 
    SeqIdxType GetOverlapPosition(const int i) const    {
        assert(i >= 0 && i <= 3);
        return ov_pos_[i];
    }

    // prints the edge information 
    void PrintInfo(void)    {
        std::cout << "Printing GraphEdgeType object info..." << std::endl;
        std::cout << "overlap length: " << overlap_ << std::endl;
        std::cout << "coverage: " << cov_ << std::endl;
        std::cout << "is visited: " << visited_ << std::endl;
        std::cout << "is resolved: " << resolved_ << std::endl;
        std::cout << "is reverse complement: " << is_rc_ << std::endl;
        return;
    }

    friend class GraphNodeType;
    friend class GraphEssential;
    friend class GraphPrune;
    friend class GraphTraversal;
    friend class AssemblyGraph;

  protected:
    SeqIdxType overlap_;        // the overlap length between the adjacent nodes
    SeqIdxType ov_pos_[4];      // the positions of the overlap (source begin, source end, target begin, target end)
    std::string cigar_;         // the cigar representation of the alignment of the overlapped region
    int num_mismatch_;          // the number of mismatches and indels found in the overlapped region
    CoverageType cov_;          // the coverage of the edge
                                // TODO: need to define the coverage of an edge using the coverages of the end nodes
    bool visited_;              // the tag indicating whether the node has been visited
    bool resolved_;             // the tag indicating whether the node has been verified (e.g., with oriention determined)
    bool src_orientation_;      // the orientation of the source node; true for positive strand, false for negative strand 
    bool is_rc_;                // the tag indicating whether the adjacent reads overlap with reverse complemetarity; true for yes, false for no 
};


#endif      // __GRAPHEDGETYPE_H_
