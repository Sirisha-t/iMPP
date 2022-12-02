#!/usr/bin/python
import sys

if len(sys.argv) == 2:
	fdir = sys.argv[1]
else:
	print("call by\n    rename.py <DBG>\n")
	exit()

name = fdir[0:fdir.rfind('.')]
fout = open(name+'.namemap.txt', 'w')
outfile = open(name+'.re.60.faa', 'w')

fin = open(fdir, 'r')
line = fin.readline()
sequences = []
headers = []
Edge_id = {}
seq = ""
i = 0
while line:
	if line[0] == '>':
		cid = line[1:]
		cid = line.rstrip()
		seq = fin.readline()
		seq = seq.rstrip()
		if(len(seq) >= 60):
			outfile.write('>PL_'+str(i))
			outfile.write('\n')
			outfile.write(seq)
			outfile.write('\n')
			fout.write(cid+'\tPL_'+str(i)+'\n')
			i += 1
	line = fin.readline()

fin.close()
fout.close()
outfile.close()
