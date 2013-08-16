import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
from operator import itemgetter
import time


#expects following format: EntrezID \t chrom \t strand \t RUD-basic \t RUDs seperated by commas 
def compareRUDs (ActiveRUDs, RestingRUDs, outfile):
	active = open(ActiveRUDs, "r")
	resting = open(RestingRUDs, "r")
	output = open(outfile, "w")
	
	outline = "# EntrezID RUD_basic_ratio RUDs_ratios [-2 = ineligible due to insignificant data] [-1 = infinity]\n"
	output.write(outline)
	
	#this is to skip the first line
	A = active.readline()
	R = resting.readline()
		
	A = active.readline()
	R = resting.readline()	
	while A != '' and R != '':
		a = A.strip('\n')
		r = R.strip('\n')
		RUD_list_A = a.split('\t')
		RUD_list_R = r.split('\t')
		if RUD_list_A[0] == RUD_list_R[0]:
			RUDs_A = RUD_list_A[4].split(',')
			RUDs_R = RUD_list_R[4].split(',')
		
			n = float(RUD_list_A[3])
			d = float(RUD_list_R[3])
		
			if n < 0 or d < 0:
				outline = RUD_list_A[0] + '\t' + "-2" + '\t' + "-2" + "\n"
			else:
		
				if d == 0.0:
					RUD_basic_ratio = "-1"
				else:
					RUD_basic_ratio = str(n/d)
		
		
				RUDs_ratio = []
				for i in range(len(RUDs_A)):
			
					n = float(RUDs_A[i])
					d = float(RUDs_R[i])			
			
					if d == 0.0:
						RUDs_ratio.append("-1")
					else:
						RUDs_ratio.append(str(n/d))
			
				outline = RUD_list_A[0] + '\t' + RUD_basic_ratio + '\t' + ",".join(RUDs_ratio) + "\n"
			
			output.write(outline)
			
			A = active.readline()
			R = resting.readline()
		else:
			print "Both input files must have the same length and use the following format:\n EntrezID RUD_basic comma_seperated_RUDs\n"
	

def main (argv):
	parser = OptionParser()
	parser.add_option("-a", "--activeRUDs", action="store", type="string", dest="ActiveRUDs", help="input Calculate3UTRsUsingStrandSpecificRNASeq RUD results for active cells", metavar="<file>")
	parser.add_option("-r", "--restingRUDs", action="store", type="string", dest="RestingRUDs", help="input Calculate3UTRsUsingStrandSpecificRNASeq RUD results for resting cells", metavar="<file>")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="outfile name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
	
	startTime = time.time()	
	
	compareRUDs(opt.ActiveRUDs, opt.RestingRUDs, opt.outfile)

	print "it took", time.time() - startTime, "seconds."

if __name__ == "__main__":
	main(sys.argv)