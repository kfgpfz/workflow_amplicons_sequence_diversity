import sys, os, csv, re
import pysam
inlines = []
A = 0
T = 0
G = 0
C = 0



samfile = pysam.Samfile(sys.argv[3], "rb")
for pileupcolumn in samfile.pileup( 'NC_008463.1', int(sys.argv[1])-1, int(sys.argv[1]),max_depth=100000):
	for pileupread in pileupcolumn.pileups:
		if pileupcolumn.pos == int(sys.argv[1])-1:
			if pileupread.query_position != None:
				#print(pileupread.query_position)
				base = pileupread.alignment.query_sequence[pileupread.query_position]
		
				if base == 'A':
					A = A + 1
				elif base == 'T':
					T = T + 1
				elif base == 'C':
					C = C + 1
				elif base == 'G':
					G = G + 1	


inlines.append([str(sys.argv[2])+":"+str(sys.argv[1])])
inlines.append(["A"]+[str(A)])
inlines.append(["T"]+[str(T)])
inlines.append(["C"]+[str(C)])
inlines.append(["G"]+[str(G)])
inlines.append(['------------'])

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "Results.txt"),"a") as outfile:
	writer = csv.writer(outfile, delimiter='	', quoting=csv.QUOTE_MINIMAL)
	for ausgabe in inlines:
		writer.writerow(ausgabe)

