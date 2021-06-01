#!/usr/bin/python

CLASSES=[[5,10,16], [15,21,26], [7,18], [11,17], [13,24], [14,20]]

import sys

def show_help():
	print('''
This script generates rules.txt file given a partial collision matrix.
The rest of the matrix is defined according to the normalization rule.
Each all undefined entities in a class are defined equal.

Usage:
	''', sys.argv[0], '''<state1>-<state2>=<prob> <state3>-<state4>=<prob> ...

Example:
	''', sys.argv[0], '''5-10=0.3 10-5=0.3 16-5=0.2
''')

# TODO
# warn if semi-detailed not fullfilled
# err if out-of-class non-zero defined
# check state is valid

if len(sys.argv)==1 or '--help' in sys.argv:
	show_help()
	sys.exit()

matrix=dict()

for rule in sys.argv[1:]:
	pair, prob=rule.split('=')
	s1, s2=list(map(int, pair.split('-')))
	if s1 not in matrix:
		matrix[s1]=dict()
	try:
		matrix[s1][s2]=float(prob)
	except ValueError:
		raise Exception('Illegal rule format', rule)
	assert 0<=float(prob)<=1

for cls in CLASSES:
	for s1 in cls:
		if s1 not in matrix:
			matrix[s1]=dict()
			for s in cls[:-1]:
				matrix[s1][s]=1.0/len(cls)
		sum=0
		for s2 in cls:
			if s2 in matrix[s1]:
				sum+=matrix[s1][s2]
		if sum>1:
			raise Exception('Normalization rule violated', sum, cls)
		
		if len(cls)>len(matrix[s1]):
			others=(1-sum)/(len(cls)-len(matrix[s1]))
			for s2 in cls:
				if s2 not in matrix[s1]:
					matrix[s1][s2]=others

for cls in CLASSES:
	print('\t'+'\t'.join(map(str, cls)))
	for s1 in cls:
		print(s1, end=' ')
		for s2 in cls:
			print('\t', str(matrix[s1][s2]), end=' ')
		print()
	print()
