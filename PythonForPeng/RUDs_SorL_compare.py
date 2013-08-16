import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
from operator import itemgetter
import time


#expects following format: EntrezID \t chrom \t strand \t RUD-basic \t RUDs seperated by commas [no spaces] 
def compareSorLfiles (file1, file2, outfile):
	f1 = open(file1, "r")
	f2 = open(file2, "r")
	output = open(outfile, "w")
	
	outline = "# EntrezID [s = shortened, l = lengthened, n = neither, i = invalid data; first letter from: " + file1 +  ", second letter from: " + file2 + "]\n"
	output.write(outline)
	
	#this is to skip the first line
	A = f1.readline()
	B = f2.readline()
		
	A = f1.readline()
	B = f2.readline()	
	while A != '' and B != '':
		a = A.strip('\n')
		b = B.strip('\n')
		RUD_list_A = a.split('\t')
		RUD_list_B = b.split('\t')
		if RUD_list_A[0] == RUD_list_B[0]:
			
			# output is the combined information from the two input files
			outline = RUD_list_A[0] + '\t' + RUD_list_A[1] + RUD_list_B[1] + "\n"
			
			output.write(outline)
			
			A = f1.readline()
			B = f2.readline()
		else:
			print "Both input files must have the same length and use the following format:\n EntrezID\tcharacter\n"
	

def main (argv):
	parser = OptionParser()
	parser.add_option("-i", "--SorL1", action="store", type="string", dest="SorL1", help="input RUDs_SorL results from one trial", metavar="<file>")
	parser.add_option("-k", "--SorL2", action="store", type="string", dest="SorL2", help="input RUDs_SorL results from a second trial", metavar="<file>")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="outfile name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
	
	startTime = time.time()	
	
	compareSorLfiles(opt.SorL1, opt.SorL2, opt.outfile)

	print "it took", time.time() - startTime, "seconds."

if __name__ == "__main__":
	main(sys.argv)