#!/usr/bin/python
import math
import re
import sys
import string
import time


def InputFileReader( contig_in, dmd_in,  read_in):
	fin = open(contig_in, 'r')
	outd = read_in[0:read_in.rfind('.')]
	outfile = open (outd+'.mapped.fasta', 'w')
	line = fin.readline()
	contig60 = set()
	mappedRead = set()
	count = 0

	##reading assembled contig fasta file
	while line:
		if line[0] == '>':
			sequence = ""
			cid = line[4:]
			cid = cid.rstrip()
		else:
			seq = line
			seq = seq.rstrip()
			if(len(seq) > 60):
				contig60.add(cid)
		line = fin.readline()

	##reading diamond blastx align file
	din = open(dmd_in, 'r')
	for line in din.readlines():
		line = line.rstrip()
		fields = line.split('\t')
		rname = fields[0]
		cname = fields[1]
		cname = cname[3:]
		pident = fields[2]
		align_len = fields[3]

		if cname in contig60:
			if float(pident) > 90.00 and int(align_len) > 20:
				mappedRead.add(rname)
				#print rname
	print len(mappedRead)
	##reading translated reads fasta file
	rin = open(read_in,'r')
	line = rin.readline()
	while line:
		if line[0] == '>':
			sequence = ""
			rid = line[1:]
			rid = rid.rstrip()
		else:
			seq = line
			seq = seq.rstrip()
			if rid in mappedRead:
				outfile.write('>'+rid)
				outfile.write('\n')
				outfile.write(seq)
				outfile.write('\n')
		line = rin.readline()


	outfile.close()

if __name__ == "__main__":
	if len(sys.argv) == 4:
		contig_fasta = sys.argv[1]
		plass_dmd = sys.argv[2]
		unpred_fasta = sys.argv[3]
		InputFileReader( contig_fasta, plass_dmd, unpred_fasta)
