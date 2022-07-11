import math
import re
import sys
import string
import time
from Bio.Seq import Seq


## DNA codon table for bacteria and archae ( genetic code : 11)
aa_dict = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S",
			  "TCA":"S","TCG":"S", "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
			  "TGT":"C","TGC":"C","TGA":"*","TGG":"W", "CTT":"L","CTC":"L",
			  "CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
			  "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R",
			  "CGA":"R","CGG":"R", "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
			  "ACT":"T","ACC":"T","ACA":"T","ACG":"T", "AAT":"N","AAC":"N",
			  "AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
			  "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A",
			  "GCA":"A","GCG":"A", "GAT":"D","GAC":"D","GAA":"E",
			  "GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}


#global out_sixfr
#out_sixfr = ""

def translate(offset, seq, out):
	out_sixfr = ""
	for x in range(int(offset), len(seq), 3):
		codon = seq[x:x + 3]
		if len(codon) < 3:
			break
		elif "N" in codon:
			#out.write("X")
			out_sixfr += 'X'
		else:
			#out.write("%s" % aa_dict[codon])
			out_sixfr += aa_dict[codon]
	return out_sixfr


def sixFrameTranslate(seqs_to_be_translated, allpredset, readinfo, out, predout):
	seq_dict = {}
	seqs_to_be_translated = seqs_to_be_translated.rstrip()
	lines = seqs_to_be_translated.split('\n')
	out_sixfr = ""
	count = 0
	for line in lines:
		count+=1
		#print count
		if line.startswith(">"):
			header = line.strip()
			seq_dict[header] = ""
		else:
			seq_dict[header] += line.strip()
	# Goes through each sequence and calls 'translate' for all reading frames
	for item in seq_dict:
		seq = seq_dict[item].upper()
		seq_rc = Seq(seq_dict[item].upper())
		seq_rc = seq_rc.reverse_complement()
		seq_rc = str(seq_rc)
		out_sixfr += item+'_frame1\n';
		#out.write("%s_frame1\n" % item)
		out_sixfr += translate(0, seq, out)
		#out.write("\n")
		out_sixfr += '\n';
		out_sixfr += item+'_frame2\n';
		#out.write("%s_frame2\n" % item)
		out_sixfr += translate(1, seq, out)
		#out.write("\n")
		out_sixfr += '\n';
		out_sixfr += item+'_frame3\n';
		#out.write("%s_frame3\n" % item)
		out_sixfr += translate(2, seq, out)
		#out.write("\n")
		out_sixfr += '\n';
		out_sixfr += item+'_frame4\n';
		#out.write("%s_frame4\n" % item)
		out_sixfr += translate(0, seq_rc, out)
		#out.write("\n")
		out_sixfr += '\n';
		out_sixfr += item+'_frame5\n';
		#out.write("%s_frame5\n" % item)
		out_sixfr += translate(1, seq_rc, out)
		#out.write("\n")
		out_sixfr += '\n';
		out_sixfr += item+'_frame6\n';
		#out.write("%s_frame6\n" % item)
		out_sixfr += translate(2, seq_rc, out)
		out_sixfr += '\n';
		#out.write("\n")

	#print out_sixfr
	## Checking if any one of the frames is longer than 20 amino acids
	if out_sixfr != "":
		header_info = set()
		lines = out_sixfr.split('\n')
		for line in lines:
			if line.startswith(">"):
				read_id = line.strip()
				read_id = read_id[1:]
				header = read_id[0:read_id.find("_")]
				#aaseq_dict[header] = ""
			else:
				seq = line.strip()
				if seq.find('*'):
					prefix = seq[0:seq.find('*')]
					if len(prefix) > 20:
						if header not in allpredset and header not in header_info :
							header_info.add(header)
							out.write(">"+header+'\n'+readinfo[header]+'\n')
							predout.write(">"+header+'\n'+readinfo[header]+'\n')
				else:
					if len(seq) > 20:
						if header not in allpredset and header not in header_info :
							header_info.add(header)
							out.write(">"+header+'\n'+readinfo[header]+'\n')
							predout.write(">"+header+'\n'+readinfo[header]+'\n')


def getSequenceInfo(read_Predictions, EP_predictions, input_read, outd):
	global seqs_to_be_translated
	global gl_longest_orf
	fin = open(input_read, 'r')
	outfiletag = input_read[0:input_read.rfind('.')]
	translatedReadsOutFile = open (outd+'/reads.translated.fasta', 'w')
	predictedReadsOutfile = open (outd+'/reads.filter.fasta', 'w')
	sequence = ""
	seqs_to_be_translated = ""
	predCount = 0
	unpredCount = 0
	allPredictedReadSet = set()
	allPredictedReadSet = read_Predictions.union(EP_predictions)
	print len(allPredictedReadSet)
	line = fin.readline()
	ReadInfo = {}
	while line:
		if line[0] == '>':
			sequence = ""
			read = line[1:]
			read = read.rstrip()
			#read_tag = read[0:read.rfind('/')]
			if(read in allPredictedReadSet):
				predCount+=1
				seq = fin.readline()
				seq = seq.rstrip()
				ReadInfo[read] = seq
				out_header = ">"
				out_header += read
				#sequence = out_header + "\n" + seq
				predictedReadsOutfile.write(out_header)
				predictedReadsOutfile.write('\n')
				predictedReadsOutfile.write(seq)
				predictedReadsOutfile.write('\n')
			else:
				unpredCount+=1
				seq = fin.readline()
				seq = seq.rstrip()
				ReadInfo[read] = seq
				out_header = ">"
				out_header += read
				sequence = out_header + '\n' + seq + '\n'

			seqs_to_be_translated += sequence
		line = fin.readline()

	if seqs_to_be_translated != "" :
		print "predicted :"+str(predCount)
		print "unpredicted :"+str(unpredCount)

		sixFrameTranslate(seqs_to_be_translated, allPredictedReadSet, ReadInfo, translatedReadsOutFile, predictedReadsOutfile)

def InputFileReader(read_gff, edge_gff, path_gff,  bwa_edgename, bwa_pathname, input_read, out_name):
	# Reading FGS predictions on input reads
	gff1 = open(read_gff, 'r')
	gff2 = open(edge_gff, 'r')
	gff3 = open(path_gff, 'r')
	readNameMap = {}
	read_Predictions = set()
	edgePred = set()
	pathPred = set()
	## reading reads gff file
	for line in gff1.readlines():
		if line[0] != '#':
			line = line.rstrip()
			fields = line.split("\t")
			rname = fields[0];
			#rname = rname[0:rname.rfind('/')]
			read_Predictions.add(rname)
	#print " Completed reading reads gff file \n"
	#print len(read_Predictions)
	## Reading egde gff file
	for line in gff2.readlines():
		if line[0] != '#':
			line = line.rstrip()
			fields = line.split("\t")
			ename = fields[0];
			edgePred.add(ename)
	#print " Completed reading edge gff file \n"
	#print len(edgePred)
	## reading paths gff file
	for line in gff3.readlines():
		if line[0] != '#':
			line = line.rstrip()
			fields = line.split("\t")
			pname = fields[0];
			pathPred.add(pname)
	#print " Completed reading paths gff file \n"
	#print len(pathPred)
	# Reading bwa mapping of reads to edges file
	bwaedge = open(bwa_edgename,'r')
	#print bwaedge
	EP_predictions = set()
	for line in bwaedge.readlines():
		if(line[0] != '@' and line[1] != 'S'):
			#print line
			line = line.rstrip()
			fields = line.split("\t")
			edge_read = fields[0];
			flag = fields[2]
			#print flag
			if flag != '*' and flag in edgePred:
				EP_predictions.add(edge_read)
	#print " Completed reading bwa edge file \n"
	#print len(EP_predictions)
	# Reading bwa mapping of reads to paths file
	bwapath = open(bwa_pathname,'r')
	for line in bwapath.readlines():
		if(line[0] != '@' and line[1] != 'S'):
			line = line.rstrip()
			fields = line.split("\t")
			path_read = fields[0];
			flag = fields[2]
			if flag != '*' and flag in pathPred:
				  EP_predictions.add(path_read)
	#print " Completed reading bwa path file \n"
	#print len(EP_predictions)
	if len(read_Predictions) != 0 and len(EP_predictions) != 0:
		print "Completed reading input file. Translating the reads now..\n"
		getSequenceInfo(read_Predictions, EP_predictions, fasta_read, out_name)

	gff1.close()
	gff2.close()
	gff3.close()
	bwaedge.close()
	bwapath.close()



if __name__ == "__main__":
	if len(sys.argv) == 8:
		read_gff_name  = sys.argv[1]
		edge_gff_name = sys.argv[2]
		path_gff_name = sys.argv[3]
		bwa_edgename = sys.argv[4]
		bwa_pathname = sys.argv[5]
		fasta_read = sys.argv[6]
		out_dirname = sys.argv[7]

		InputFileReader(read_gff_name, edge_gff_name, path_gff_name, bwa_edgename, bwa_pathname, fasta_read, out_dirname)
