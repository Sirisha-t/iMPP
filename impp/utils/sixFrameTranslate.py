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
		elif codon in aa_dict:
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
		if line.startswith(">"):
			header = line.strip()
		else:
			seq = line.strip()
	# For each sequence call translate on all reading frames
	seq = seq.upper()
	seq_rc = Seq(seq)
	seq_rc = seq_rc.reverse_complement()
	seq_rc = str(seq_rc)
	out_sixfr += header+'_frame1\n';
	out_sixfr += translate(0, seq, out)
	out_sixfr += '\n';
	out_sixfr += header+'_frame2\n';
	out_sixfr += translate(1, seq, out)
	out_sixfr += '\n';
	out_sixfr += header+'_frame3\n';
	out_sixfr += translate(2, seq, out)
	out_sixfr += '\n';
	out_sixfr += header+'_frame4\n';
	out_sixfr += translate(0, seq_rc, out)
	out_sixfr += '\n';
	out_sixfr += header+'_frame5\n';
	out_sixfr += translate(1, seq_rc, out)
	out_sixfr += '\n';
	out_sixfr += header+'_frame6\n';
	out_sixfr += translate(2, seq_rc, out)
	out_sixfr += '\n';

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
						if header not in header_info :
							header_info.add(header)
							out.write(">"+header+'\n'+readinfo[header]+'\n')
							predout.write(">"+header+'\n'+readinfo[header]+'\n')
				else:
					if len(seq) > 20:
						if header not in header_info :
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
	count = 0
	allPredictedReadSet = set()
	allInputReadSet = set()
	allPredictedReadSet = read_Predictions.union(EP_predictions)
	line = fin.readline()
	ReadInfo = {}
	start=time.time()
	while line:
		if line[0] == '>':
			sequence = ""
			read = line[1:]
			read = read.rstrip()
			seq = fin.readline()
			seq = seq.rstrip()
			ReadInfo[read] = seq
			if read in allPredictedReadSet:
				out_header = ">"
				out_header += read
				predictedReadsOutfile.write(out_header)
				predictedReadsOutfile.write('\n')
				predictedReadsOutfile.write(seq)
				predictedReadsOutfile.write('\n')

			allInputReadSet.add(read)

		line = fin.readline()

	UnPredReadSet = allInputReadSet.difference(allPredictedReadSet)

	end=time.time()
	#print("reading fasta to get unpred set in + getting filtered reads: ")
	#print((end-start))

	start=time.time()
	for read in UnPredReadSet:
		count = count + 1
		out_header = ">"
		out_header += read
		sequence = out_header + '\n' + ReadInfo[read]
		sixFrameTranslate(sequence, allPredictedReadSet, ReadInfo, translatedReadsOutFile, predictedReadsOutfile)
		#seqs_to_be_translated += sequence
	end=time.time()



def InputFileReader(read_gff, edge_gff, path_gff,  bwa_edgename, bwa_pathname, input_read, out_name):
	# Reading FGS predictions on input reads
	gff1 = open(read_gff, 'r')
	gff2 = open(edge_gff, 'r')
	gff3 = open(path_gff, 'r')
	readNameMap = {}
	read_Predictions = set()
	edgePred = {}
	pathPred = {}
	## reading reads gff file
	start = time.time()
	for line in gff1.readlines():
		if line[0] != '#':
			line = line.rstrip().split('\t')
			rname = line[0];
			#rname = rname[0:rname.rfind('/')]
			read_Predictions.add(rname)
	end = time.time()
	#print(" Completed reading reads gff file in :")
	#print((end-start))
	#print(len(read_Predictions))

	## Reading egde gff file
	start = time.time()
	for line in gff2.readlines():
		if line[0] != '#':
			line = line.rstrip().split('\t')
			ename = line[0];
			fields = ename.rstrip().split(',')
			ename = fields[0];
			#print ename
			edgePred[ename] = 1;
	end = time.time()
	#print(" Completed reading edges gff file in :")
	#print((end-start))
	print ("Completed reading edge gff file \n")
	#print len(edgePred)
	
	## reading paths gff file
	start = time.time()
	for line in gff3.readlines():
		if line[0] != '#':
			line = line.rstrip().split('\t')
			pname = line[0];
			fields = pname.rstrip().split(',')
			pname = fields[0];
			pathPred[pname] = 1;
	end = time.time()
	#print(" Completed reading paths gff file in :")
	#print((end-start))
	print (" Completed reading paths gff file \n")
	#print len(pathPred)
	
	# Reading bwa mapping of reads to edges file
	bwaedge = open(bwa_edgename,'r')
	#print bwaedge
	EP_predictions = set()
	start = time.time()
	for line in bwaedge.readlines():
		if(line[0] != '@' and line[1] != 'S'):
			line = line.rstrip().split('\t')
			edge_read = line[0];
			flag = line[2]
			if flag != '*':
				flag = flag.rstrip().split(',')
				ename = flag[0];
				if ename in edgePred:
					EP_predictions.add(edge_read)
	end = time.time()
	#print(" Completed reading bwa-edge file in :")
	#print((end-start))
	print (" Completed reading bwa edge file \n")
	#print len(EP_predictions)
	
	# Reading bwa mapping of reads to paths file
	start = time.time()
	bwapath = open(bwa_pathname,'r')
	for line in bwapath.readlines():
		if(line[0] != '@' and line[1] != 'S'):
			line = line.rstrip().split('\t')
			path_read = line[0];
			flag = line[2]
			if flag != '*':
				flag = flag.rstrip().split(',')
				pname = flag[0];
				if pname in pathPred:
					EP_predictions.add(edge_read)
	end = time.time()
	#print(" Completed reading bwa-path file in :")
	#print((end-start))
	print (" Completed reading bwa path file \n")
	#print len(EP_predictions)
	
	if len(read_Predictions) != 0 and len(EP_predictions) != 0:
		#print("Completed reading input file. Translating the reads now..\n")
		start = time.time()
		getSequenceInfo(read_Predictions, EP_predictions, fasta_read, out_name)
		end = time.time()
		#print(" Completed running getSequenceInfo function in :")
		#print((end-start))


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
