#!/usr/bin/python
import re
import sys

if __name__ == "__main__":
    input_file = sys.argv[1]
    freq = sys.argv[2]
    print freq
    #out_dirname = sys.argv[3]
    global seqs_to_be_translated
    global gl_longest_orf
    fin = open(input_file, 'r')
    outfiletag = input_file[0:input_file.rfind('.')]
    shortoutfile = open (outfiletag+'.short.'+freq+'.fq', 'w')
    longoutfile = open (outfiletag+'.long.'+freq+'.fq', 'w')
    sequence = ""
    line = fin.readline()
    while line:
        if line[0] == '>':
            sequence = ""
            read = line[0:]
            read = read.rstrip()
            #print read
        else:
                #seq = fin.readline()
                seq = line
                seq = seq.rstrip()
                #print len(seq)

                if(len(seq) > int(freq)):
                    #print len(seq)
                    longoutfile.write(read)
                    longoutfile.write('\n')
                    longoutfile.write(seq)
                    longoutfile.write('\n')
                else:
                    shortoutfile.write(read)
                    shortoutfile.write('\n')
                    shortoutfile.write(seq)
                    shortoutfile.write('\n')
        line = fin.readline()
