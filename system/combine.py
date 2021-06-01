#!/usr/bin/python

# Scans '$1' forder for *.xls files and puts them in alphabetical order into one file '$2' line by line

import sys, os

xls_dir=sys.argv[1]
out_path=sys.argv[2]

files=[f for f in os.listdir(xls_dir) if f.endswith('.xls')]

fl=open(out_path, 'w')

for f in files:
	lines=open(os.path.join(xls_dir, f)).readlines()
	for l in lines:
		j, val=l.split('\t')
		fl.write(val.strip() + '\t')
	fl.write('\n')

fl.close()
