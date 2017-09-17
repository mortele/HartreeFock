import sys


inFileName  = "Integrator/IntegralTables/coulomb.dat";
outFileName = "Integrator/IntegralTables/coulomb2.dat";

with open(inFileName, 'r') as inFile :

	outFile = open(outFileName, 'w')

	for line in inFile :
		l = line.split()
		outFile.write("%d %d %d %d %f\n" % (int(l[0])-1, int(l[1])-1, int(l[2])-1, int(l[3])-1, float(l[4])))

	outFile.close()

