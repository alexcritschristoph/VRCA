import sys
import json
import os
with open(sys.argv[2]) as data_file:
		data = json.load(data_file)

with open(sys.argv[1]) as virus_host:
		for line in virus_host:
			try:
				name = line.split("\t")[4].split()[0]
				for key in data.keys():
					if name in key.split():
						print line.split("\t")[0].strip() + "," + line.split("\t")[4].strip() + "," + key
						os.system("wget 'http://www.ebi.ac.uk/ena/data/view/" + line.split("\t")[0].strip()+"&display=fasta'")
						break
			except:
				print "err"
