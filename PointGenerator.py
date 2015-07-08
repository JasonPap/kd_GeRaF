#!/usr/bin/python

import random
import sys

# Get arguments
size = int(sys.argv[1])
dimension = int(sys.argv[2])
outputFile = sys.argv[3]

# Open output file
out = open(outputFile, "w")

# Write random Points to file
for line in range(0, size):
	for coordinate in range(0, dimension):
		out.write("%s " % random.uniform(0,100))
	out.write("\n")

out.close()


