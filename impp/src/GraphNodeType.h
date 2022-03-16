#ifndef __GRAPHNODETYPE_H_
#define __GRAPHNODETYPE_H_

#include <iostream>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <assert.h>

#include "DataType.h"   

// The assembly graph node class
class GraphNodeType {
  public:
    
    // empty constructor function
    GraphNodeType(void) {
        str_ = nullptr;
        str_len_ = 0;
        cov_ = 0;
        visited_ = resolved_ = false;
        orientation_ = true;
    }

    // constructor function from a char array
    explicit GraphNodeType(const char *s)  {
        str_len_ = strlen(s);
        str_ = new char [str_len_ + 1];
        strcpy(str_, s);
        cov_ = 0;
        visited_ = resolved_ = false;
        orientation_ = true;
        return;
    }

    // destructor
    ~GraphNodeType(void)    {
        if(str_len_ > 0)    delete []  str_;
    }

    // assignment operator
    GraphNodeType& operator=(const GraphNodeType &n) {
        this->str_len_ = n.GetSeqLen();
        this->str_ = new char [this->str_len_ + 1];
        strcpy(this->str_, n.str_);
        this->cov_ = n.cov_;
        this->visited_ = n.visited_;
        this->resolved_ = n.resolved_;
        this->orientation_ = n.orientation_;
        return *this;
    } 

    // updating the sequence of a node
    void SetSequence(char *s)    {
        str_len_ = strlen(s);
        str_ = new char [str_len_ + 1];
        strcpy(str_, s);
        return;
    }

    // access the length of the sequence contained in the node
    SeqIdxType GetSeqLen(void) const {
        return str_len_;
    } 

    // TODO: to be updated with BioSequence class data access
    const char* GetStrPtr(void) const {
        return str_;
    }

    // setting the coverage of the current node
    void SetCoverage(const CoverageType c)    {
        cov_ = c;
        return;
    }

    // returns the coverage of the current node
    CoverageType GetCoverage(void) const    {
        return cov_;
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

    // returns the orientation of the read
    bool GetOrientation(void) const  {
        return orientation_;
    }

    // setting the node's orientation
    void SetOrientation(const bool b) {
        orientation_ = b;
        return;
    }

    // print the node information
    void PrintInfo(void)    {
        std::cout << "Printing GraphNodeType object info... " << id_ << std::endl;
        std::cout << "sequence: " << str_ << std::endl;
        std::cout << "length: " << str_len_ << std::endl;
        std::cout << "coverage: " << cov_ << std::endl;
        std::cout << "is visited: " << visited_ << std::endl;
        std::cout << "is resolved: " << resolved_ << std::endl;
        std::cout << "orientation: " << orientation_ << std::endl;
        return;
    }

    friend class GraphEdgeType;
    friend class GraphEssential;
    friend class GraphPrune;
    friend class GraphTraversal;
    friend class AssemblyGraph;

  protected:
    IDType id_;                 // DEBUG
    char *str_;                 // the char array that holds the sequence
    SeqIdxType str_len_;        // the length of the sequence
    CoverageType cov_;          // the coverage of the node
                                // TODO: to define profile type, waiting for BioSequence definition
                                // Or, can we compute but not actually store the frequency informaiton?
    bool visited_;              // the tag indicating whether the node has been visited
    bool resolved_;             // the tag indicating whether the node has been verified (e.g., with oriention determined)
    bool orientation_;          // the tag indicating the orientation of the read; true for plus strand, false for minus strand

};

#endif  // __GRAPHNODETYPE_H_
