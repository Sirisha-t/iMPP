import sys, subprocess, gzip


def convertFq(fastq, outd):
	line_n =0
	line_buffer = 0
	line_id = 1
        read_id = 0
	readMap = {}
	outfile = open(outd+'/reads.fq2fa.fasta', 'w')
	fqoutfile = open(outd+'/reads.re.fastq', 'w')
	seqNameMap = open(outd+'/reads.namemap.txt', 'w')
	fastas = 1
	fasta_length =0
	for line in fastq:
		line_n += 1
			#line_id += 1
		if line_n == 10000:
			line_buffer += 10000
			line_n =0
		if line_id == 4:
			fqoutfile.write(line)
			line_id = 1
		elif line_id == 3:
			fqoutfile.write(line)
			line_id += 1
		elif line_id == 2:
			line_id += 1
			fasta_line = line
			fasta_length += len(fasta_line.strip())
			fqoutfile.write(line)
			outfile.write(fasta_line)
			fastas += 1
		else:
			if '@' not in line:
				print "ERROR: Input file is not in correct FASTQ format.\n"
				break
			else:
				rname = line[1:]
                		readMap[rname] = read_id;
                		fasta_header = ">"
                		fasta_header += str(readMap[rname])
				line_id += 1
 		        	read_id += 1
				fastq_header = fasta_header.replace('>', '@')
				fqoutfile.write(fastq_header)
				fqoutfile.write("\n")
				outfile.write(fasta_header)
				outfile.write("\n")
	#sorted_readMap = sorted(readMap.items(), key=lambda x: x[1])
	seqNameMap.write("### Read name \t Read id ###\n")
	for name in readMap:
		outputline = name.rstrip()
		outputline += "\t" 
		outputline +=  str(readMap[name])
		#print outputline
		seqNameMap.write(outputline)
		seqNameMap.write("\n")	
	outfile.close()
	#print 'FASTA records written', fastas, 'average length of fasta sequences ', float(fasta_length//fastas)

if __name__ == '__main__':
	infile = sys.argv[1]
	outd = sys.argv[2]
	fasta = infile[0:infile.rfind('.')]
	fastq = open(infile, 'r')

	convertFq(fastq, outd)
