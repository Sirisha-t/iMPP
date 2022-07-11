#!/usr/bin/python
import math
import re
import sys
import string
import time


def InputFileReader(filt_in, trans_in):
	fin = open(trans_in, 'r')
	outd = filt_in[0:filt_in.rfind('.')]
	outfile = open (outd+'.pred.fasta', 'w')
	#line = fin.readline()
	count = 0
	trans = set()
	##reading trnaslated fasta file
	for line in fin.readlines():
		if line[0] == '>':
			sequence = ""
			cid = line[1:]
			cid = cid.rstrip()
			trans.add(cid)

	filtin = open(filt_in, 'r')
	line = filtin.readline()
	while line:
		if line[0] == '>':
			sequence = ""
			cid = line[1:]
			cid = cid.rstrip()
			if cid not in trans:
				seq = filtin.readline()
				seq = seq.rstrip()
				outfile.write('>'+str(cid))
				outfile.write('\n')
				outfile.write(seq)
				outfile.write('\n')
		line = filtin.readline()
		
	outfile.close()

if __name__ == "__main__":
	if len(sys.argv) == 3:
		filter_fasta = sys.argv[1]
		trans_fasta = sys.argv[2]

		InputFileReader(filter_fasta, trans_fasta)
