#!/usr/bin/python

# USER CONFIG

TITLE="0.45-0.45-0.1 (4R)"

GRID_X_PERIOD=250
GRID_Y_PERIOD=0.1

LABEL_X="Density"
LABEL_Y="Space coordinate"

# SCRIPT, DO NOT MODIFY BELOW THIS LINE

import sys, os

def plot(csv_path):
	#out_path=

files=[f for f in os.listdir() if f.endswith('.csv')]

for f in files:
	plot(f)
