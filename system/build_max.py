#!/usr/bin/python

import sys

if len(sys.argv)<2:
	print('usage:', sys.argv[0], '<file1> <file2> ... <fileN>')
	sys.exit()

for f in sys.argv[1:]:
	lines=open(f).readlines()

	top=None

	for ln in lines:
		key, val=ln.split('\t')
		val=float(val)
		if top is None or val>top:
			top=val
	print(top)
