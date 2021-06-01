#!/usr/bin/python

import sys
from functools import reduce

# ALGORITHM:
# obtain file name from sys.argv[1]
# check that file exists
# try read contents from file
# parse contents:
	# skip empty lines and comments
	# if non-empty line, read states from it
	# check states are valid
	# read next N lines
		# check line syntax:
			# has N+1 words
			# fist word is in states[]
			# other words are floats
		# add rules to table
	# check table
		# give error if normalization fails
		# give warning if semi-detailed balance fails
# generate collision function to argv[2]

def error(message):
	print()
	print("ERROR", message)
	print()
	sys.exit(1)

def extract_words(line):
	return [x for x in line.replace('\t', ' ').split(' ') if x.strip()]

def validate_prob(word, line):
	try:
		num=float(word)
		if num>=0 and num<=1:
			return
	except ValueError:
		pass
	error("Not a probability (" + str(word) + ') in line: ' + line)

def validate_state(word, line):
	try:
		num=int(word)
		if num>=0 and num<=255:
			return
	except ValueError:
		pass
	error("Not a state (" + str(word) + ') in line: ' + line)

rules={}

def init_rules():
	global rules
	for i in range(256):
		rules[i]={i: 1.0}

def add_rule(s0, s1, p):
	global rules
	rules[s0][s1]=p
	if p==0:
		del rules[s0][s1]

def is_deterministic():
	for r in rules:
		if len(rules[r])>1:
			return False
	return True

init_rules()

if len(sys.argv) != 3:
	print()
	print("Usage:", sys.argv[0], "<rules_file.txt>", "<output_file.c>")
	print()
	sys.exit()

try:
	lines=[x for x in open(sys.argv[1]).readlines() if x.strip()]
except IOError:
	error("ERROR: Failed to read from file " + sys.argv[1])

i=0
while i<len(lines):
	line=lines[i].strip()
	states=extract_words(line)
	if len(set(states)) != len(states):
		error("Duplicate state(s) in line: " + line) 
	for word in states:
		validate_state(word, line)
	if i+len(states) >= len(lines):
		error("Unexpected end of file, need more lines for states (" +\
			', '.join(states) + '), given in line ' + line);
	i+=1
	for j in range(len(states)):
		words=extract_words(lines[i+j])
		if len(words) != len(states)+1:
			error("Too few words in line " + lines[i+j])
		validate_state(words[0], lines[i+j])
		if words[0] != states[j]:
			error("Invalid state (" + words[0] + "), need (" +\
				states[j] + ") in line " + lines[i+j])
		for k in range(len(words)-1):
			validate_prob(words[k+1], lines[i+j])
			add_rule(int(words[0]), int(states[k]), float(words[k+1]))
	i+=len(states)

# check normalization
for s0 in range(256):
	sum_prob=0
	if len(rules[s0])>0:
		sum_prob=reduce(lambda x, y: x+y,
			[rules[s0][s1] for s1 in rules[s0]])
	if abs(sum_prob-1)>0.0000001:
		error("Normalization fails for state (" + str(s0) +\
			"), sum prob is " + str(sum_prob))

# generate code
try:
	f=open(sys.argv[2], 'w')
	f.write("char collide(char cell)\n")
	f.write("{\n")
	if not is_deterministic():
		f.write("\tdouble r=drand48();\n")
	f.write("\tswitch(cell) {\n")
	for s0 in range(256):
		if len(rules[s0])==1 and s0 in rules[s0]:
			continue # trivial rule
		f.write("\t\tcase " + str(s0) + ":\n")
		cur_prob=0
		for s1, prob in rules[s0].items():
			cur_prob+=prob
			if cur_prob<1.0:
				f.write("\t\t\tif (r<" + str(cur_prob) + ") return " +\
					str(s1) + ";\n")
			else:
				f.write("\t\t\treturn " + str(s1) + ";\n");
	f.write("\t\tdefault:\n")
	f.write("\t\t\treturn cell;\n")
	f.write("\t}\n")
	f.write("}\n")
except IOError:
	error("Failed to write output to file " + sys.argv[2])
