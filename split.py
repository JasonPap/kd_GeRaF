import sys


#read command line arguments
filename = sys.argv[1]
num = int(sys.argv[2])

print "Spliting file into " + str(num) + " parts"

outputfiles = []
for i in range(0,num):
	outputfiles.append(open(filename+str(i), 'w'))

with open(filename) as inputfile:
	i=0
	for line in inputfile:
		outputfiles[i%num].write(line)
		i+=1

for i in range(0,num):
	outputfiles[i].close()

print "done"