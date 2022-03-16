#ifndef __STRINGUTILS_H_
#define __STRINGUTILS_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/tokenizer.hpp>

#include "DataType.h"

namespace StringUtils   
{
    // Splits a single string into a set of strings using a specific delimiter
    // multiple delimiters are considered as a single one
    // parameters:
    //      s:  the input string to be decomposed
    //      d:  the string for concatenated delimiters; e.g., to use either space ' ' or tab '\t' as delimter set d = " \t"
    //      vs: the vector that holds the output decomposed strings
    void SplitByDelimiter(const std::string &s, const std::string &d, std::vector<std::string> &vs);

    // Gets the reverse complementary sequence of a nucleotide sequence
    // parameters:
    //      s: the input nucleotide string
    // Returns the reverse-complementary string
    std::string NARevComplementary(const std::string &s); 

    // Perform in-place reverse complementary conversion of a char array
    // parameters:
    //      s: the input char array
    //      l: the length of the char array
    void InplaceRevComp(char *s, const int l);
    
}

#endif  //__STRINGUTILS_H_
