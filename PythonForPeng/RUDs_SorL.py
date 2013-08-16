
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
from operator import itemgetter
import time


def Find_shortening_or_lengthening (infile, outfile, s-limit, l-limit):
	A_vs_R = open(infile, "r")
	output = open(outfile, "w")
	
	outline = "# EntrezID [s = shortened, l = lengthened, n = neither, i = invalid data]" + "\t" + "{S-Limit = " + str(s-limit) + ", L-Limit = " + str(l-limit) + "}" + "\n"
	output.write(outline)
	
	#this is to skip the first line
	K = A_vs_R.readline()
		
	K = A_vs_R.readline()
	while K != '':
		k = K.strip('\n')
		K_list = k.split('\t')
		J = float(K_list[1])
		
		
		if J == -2:
			X = 'I'
		elif J == -1 or J > l-limit:
			X = 'L'
		elif J < s-limit:
			X = 'S'
		else:
			X = 'N'
			
		outline = K_list[0] + '\t' + X + "\n"
						
		output.write(outline)
						
		K = A_vs_R.readline()	
			

def main (argv):
	parser = OptionParser()
	parser.add_option("-i", "--active_vs_resting", action="store", type="string", dest="Active_vs_Resting", help="input a RUDcompare output file", metavar="<file>")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="outfile name", metavar="<file>")
	#The choice of shortening-limit and lengthening-limit is very important
	parser.add_option("-s", "--shortening-limit", action="store", type="float", dest="s-limit",help="float value that Active_RUD/Rest_RUD must be less than to be considered shortening", metavar="<float>")
	parser.add_option("-l", "--lengthening-limit", action="store", type="float", dest="l-limit",help="float value that Active_RUD/Rest_RUD must be greater than to be considered lengthening", metavar="<float>")
	
	(opt, args) = parser.parse_args(argv)
	
	if len(argv) < 8:
		parser.print_help()
		sys.exit(1)
	
	startTime = time.time()	
	
	Find_shortening_or_lengthening(opt.Active_vs_Resting, opt.outfile, opt.s-limit, opt.l-limit)

	print "it took", time.time() - startTime, "seconds."

if __name__ == "__main__":
	main(sys.argv)