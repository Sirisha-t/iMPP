#include "StringUtils.h"

using namespace std;

// Splits a single string into a set of strings using a specific delimiter
// multiple delimiters are considered as a single one
// parameters:
//      s:  the input string to be decomposed
//      d:  the string for concatenated delimiters; e.g., to use either space ' ' or tab '\t' as delimter set d = " \t"
//      vs: the vector that holds the output decomposed strings
void StringUtils::SplitByDelimiter(const std::string &s, const std::string & d, vector<string> &vs)
{
    boost::split(vs, s, boost::is_any_of(d), boost::token_compress_on);
    return;
}

// Gets the reverse complementary sequence of a nucleotide sequence
// parameters:
//      s: the input nucleotide string
// Returns the reverse-complementary string
std::string StringUtils::NARevComplementary(const std::string &s) {
    // TODO: should be included into the BioSequence class
    string t = s;
    SeqIdxType l = s.length();
    //cout << "DEBUG: " << l << endl;
    for(size_t i = 0; i < l; ++ i)   {
        switch(s[i])    {
            case 'A':
                t[l - i - 1] = 'T';
                break;
            case 'C':
                t[l - i - 1] = 'G';
                break;
            case 'G':
                t[l - i - 1] = 'C';
                break;
            case 'T':
                t[l - i - 1] = 'A';
                break;
            default:
                cerr << "Warning: StringUtils::NARevComplementary:   Unrecongnized symbol, assuming \'N\'." << endl;
                t[l - i - 1] = 'N';
        }
        //cout << "DEBUG: " << i << " " << l - i - 1 << endl;
    }
    //cout << "DEBUG: END of for loop." << endl;
    return t;
}

// Perform in-place reverse complementary conversion of a char array
// parameters:
//      s: the input char array
//      l: the length of the char array
void StringUtils::InplaceRevComp(char *s, const int l)  {
    assert(strlen(s) == l); // make sure the lengths match
    char *t = new char[l + 1];  t[l] = '\0';
    for(size_t i = 0; i < l; ++ i)   {
        switch(s[i])    {
            case 'A':
                t[l - i - 1] = 'T';
                break;
            case 'C':
                t[l - i - 1] = 'G';
                break;
            case 'G':
                t[l - i - 1] = 'C';
                break;
            case 'T':
                t[l - i - 1] = 'A';
                break;
            default:
                cerr << "Warning: StringUtils::NARevComplementary:   Unrecongnized symbol, assuming \'N\'." << endl;
                t[l - i - 1] = 'N';
        }
        //cout << "DEBUG: " << i << " " << l - i - 1 << endl;
    }
    strcpy(s, t);
    delete [] t;
    return;
}