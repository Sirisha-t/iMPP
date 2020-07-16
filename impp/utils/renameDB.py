#!python

import sys

if len(sys.argv) == 2:
    fdir = sys.argv[1]
else:
    print("call by\n    renameDB.py <DBG>\n")
    exit()

name = fdir[0:fdir.rfind('.')]

fin = open(fdir, 'r')
sequences = []
headers = []
Edge_id = {}
seq = ""
i = 0
for line in fin.readlines():
    if line[0] == '>':
        sequences.append("")
        header = []
        if ':' in line:
            source = line[1:line.find(':')]
            Edge_id[source] = i
            targets = line[line.find(':')+1:-2]

            header.append(source)
            cont = targets.split(',')
            for it in cont:
                header.append(it)
        else:
            source = line[1:-2]
            Edge_id[source] = i
            header.append(source)

        i += 1
        headers.append(header)
    else:
        sequences[-1] += line[0:-1]

fin.close()

fout1 = open(name+'.namemap.txt', 'w')
fout2 = open(name+'.dbGraph.fq', 'w')
for idx in range(0,len(sequences)):
    this_header = headers[idx]
    source = this_header[0]
    fout1.write(source)
    fout1.write('\n')

    new_header = ">"
    new_header += str(Edge_id[source])

    if len(this_header) > 1:
        new_header += ':'
        new_header += str(Edge_id[this_header[1]])
        for j in range(2,len(this_header)):
            new_header += ','
            new_header += str(Edge_id[this_header[j]])

    fout2.write(new_header)
    fout2.write('\n')
    fout2.write(sequences[idx])
    fout2.write('\n')
fout1.close()
fout2.close()
