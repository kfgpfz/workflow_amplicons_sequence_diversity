#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os, csv

try:
	Fitnessergebnis = sys.argv[1]
except IndexError:
	sys.exit("Error: No input file given.")


inlines = []
inlines.append(["Position","A","T","C","G"])


with open (Fitnessergebnis, 'rt') as infile:
	#i = 0	
	for line in infile:
				
		if line.startswith("NC_008463.1:"):
			position = line.split("NC_008463.1:")[1][:-1]
			#print(position)
			inlines.append([]) # New Line
			inlines[-1].append(int(position))
			inlines[-1].append("0")	
			inlines[-1].append("0")
			inlines[-1].append("0")
			inlines[-1].append("0")		
		elif line.startswith("A	"):
			A = line.split("A")[1][:-1]
			inlines[-1].insert(1,int(A))
		elif line.startswith("T	"):
			T = line.split("T")[1][:-1]
			inlines[-1].insert(2,int(T))
		elif line.startswith("C	"):
			C = line.split("C")[1][:-1]
			inlines[-1].insert(3,int(C))
		elif line.startswith("G	"):
			G = line.split("G")[1][:-1]
			inlines[-1].insert(4,int(G))
		#if i == 14:
			#break 
		#i=i+1

#print(inlines)

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), Fitnessergebnis+"_output.txt"),"w") as outfile:
	writer = csv.writer(outfile, delimiter='	', quoting=csv.QUOTE_MINIMAL)
	for outline in inlines:
		writer.writerow(outline)

inlines = []

with open(Fitnessergebnis+"_output.txt", 'rt') as kurz:
	for line in kurz:
	
		teilstring=line.split()
		ausgabe=teilstring[0:5]
		inlines.insert(357, ausgabe)
		#print(ausgabe)
		
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "Results_bereinigt.txt"),"w") as outfile:
	writer = csv.writer(outfile, delimiter='	', quoting=csv.QUOTE_MINIMAL)
	for ausgabe in inlines:
		writer.writerow(ausgabe)
