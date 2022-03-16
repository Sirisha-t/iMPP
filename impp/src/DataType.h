#ifndef __DATATYPE_H_
#define __DATATYPE_H_

#include <iostream>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <assert.h>


//  This is the index for a large sequence database. 
//  We need to use unsigned 64bit because the database can be huge (e.g., a large metagenomic sequencing dataset). 
//  Note that using unsigned 32bit we can only handle ~4B letters, which is smaller than typical large databases.
//  Hence we are using 64bit here. 
typedef uint64_t DBIdxType;


//  This is the index for a single sequence. 
//  The worst-case scenario is that we handle a long chromosome, which is typically shorter than 4B letters. 
//  So, unsigned 32bit is fine. 
//  When handling cases where length of sequence is involved, such as sequence overlap, longest common prefix, alignment length etc., 
//  SeqIdxType should also be used.
typedef uint32_t SeqIdxType;  


//  This is the index for the read IDs. (Or other IDs if applicable.) 
//  Typically, an NGS dataset will contain less than 4B reads, although with many exceptions. 
//  Since we cannot handle that much data anyway, we will split the large dataset into smaller ones. 
//  Each small dataset should have <4B reads, which can be represented using unsigned 32bit.
typedef uint32_t IDType;


// This is the type for frequency (or percentage) definition.
// This type should be used when describing the frequency of a residue as integer percentage
// the ASCII code of the value should be subtracted by 33 to convert to the actual percentage
// for example, "+" has an ASCII value of 43, subtract 33 we have 10, indicating that the frequency is 10%
typedef uint8_t FrequencyType;


// This is the type for sequence coverage.
// Most cases we expect the coverage to be less than 4B. If it is larger, we set it as 4B.
typedef uint32_t CoverageType;

#endif      // __DATATYPE_H_