#!/usr/bin/env python
# Build a data structure: gene_id: [list of PA peaks]
# Calculate the characteristics for each gene_id

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect
from operator import itemgetter
try:
	import cPickle as pickle
except:
	import pickle

import time


sys.path.append('/Users/andrewleinbach/Desktop/SICER1.1/SICER/lib')
sys.path.append('/Users/andrewleinbach/Desktop/SICER1.1/SICER/extra')
sys.path.append('/Users/andrewleinbach/Desktop/SICER1.1/SICER/extra/tools/RNASeq')

import Utility_extended
import associate_tags_with_regions
import GenomeData
import get_total_tag_counts

import Entrez
#import Calculate3UTRUsageIndexFromCuratedGenes
from Entrez import KnownEntrezGenes
from Entrez import EntrezGene

plus = re.compile("\+")
minus = re.compile("\-")
comment = re.compile("#|track")

def Calculate3UTRUsage(entrez_genes, bedfile, chroms, outfile, threshold, PAfile, extension, index):
	"""
	entrez genes are made sure to be on one strand, 
	the bed file are reads for that strand

	entrez_genes is a KnownEntrezGenes class object
	The raw read file needs to conform to bed format

	column_index: column in bed file for sorting

	"""
	# Separate reads by chrom 
	rawreadslibName1 = (bedfile).split('/')[-1]
	rawreadssuffix1 = rawreadslibName1.split('.')[-1] 
	rawreadslibName1 = rawreadslibName1.split('.')[0]
	rawreadsextension1 = "-" + rawreadslibName1 +'.' + rawreadssuffix1 + "1"
	if Utility_extended.fileExists(bedfile):
		if Utility_extended.chrom_files_exist(chroms, rawreadsextension1) != 1:
			# Separate by chrom and sort by start
			print chroms, rawreadsextension1, " files do not exist, separate by chroms and sort each file according to the second column. "
			Utility_extended.separate_by_chrom_sort(chroms, bedfile, rawreadsextension1, str(index))
	else:
		print bedfile, " is not found"
		sys.exit(1)

	#This part is to access the polyadenylation sites
	PA1 = open(PAfile, 'r')
	
	PAsiteslist = []
	PA2 = 'i'
	while PA2 != '':
		PA2 = PA1.readline()
		if PA2 != '':
			PA3 = PA2.strip('\n')
			PA4 = PA3.split('\t')
			PAsiteslist.append((PA4[0],PA4[1]))

	PA1.close()

	# Here the output is 'a', i.e. the output is appended to an existing file instead of creating one
	outf = open(outfile, 'a')	
	for chrom in chroms: 
		if chrom in entrez_genes.chroms:
			# a KnownEntrezGenes object
			entrez_genes_by_chrom =  Entrez.KnownEntrezGenes([chrom], entrez_genes.subset_by_chrom(chrom))
			# Get the read locations
			if Utility_extended.fileExists(chrom + rawreadsextension1):
				f = open(chrom + rawreadsextension1, 'r')
				tag_positions = []
				for line in f:
					line = line.strip()
					sline = line.split()
					#make sure the extension is always 0, otherwise the rest of the program might not work as intended
					tag_positions.append(associate_tags_with_regions.tag_position(sline, 0))
				
				f.close()
				if not Utility_extended.is_list_sorted(tag_positions):
					tag_positions.sort()					
				#By this point tag_positions is a sorted list of all the reads located on the strand and chromosome the code is currently dealing with

				for entrez_id in entrez_genes_by_chrom.entrez_ids:
					gene = entrez_genes_by_chrom.entrez_genes[entrez_id] # an EntrezGene class object
					# get_3UTRs gets the ENTREZ 3'UTR, which appears to generally give the beginning of the 3'UTR and a site very close to the most distal polyadenylation site
					three_UTRs = gene.get_3UTRs()
					# Mastertuplemaker uses the ENTREZ 3'UTR and the polyA sites given to create the true data for the 3'UTR needed for CUTR_vs_AUTR to work
					true3UTRstarts, true3UTRends, UTRregion_start, UTRregion_end, UTRbeginning = Mastertuplemaker(three_UTRs,PAsiteslist,chrom,gene.strand, extension)
					#value should always be 1 as only 3'UTR with more than 1 polyA site need be considered
					if len(true3UTRends) > 1:
						#find all reads inside the 3'UTR
						inside_reads = associate_tags_with_3UTR(tag_positions, UTRregion_start, UTRregion_end)
						#finds reads in each region of the 3'UTR and calculates aUTR/cUTR for each of them
						#PolyAsites potentially useful for output
						RUDs, basic_RUD, PolyAsites = CUTR_vs_AUTR(true3UTRstarts, true3UTRends, inside_reads, gene.strand, threshold)
						
						#important if one wants to output gene_symbol information
						gene_symbol = []
						for mytranscript in gene.transcripts:
							if mytranscript.additional_annotations[0] not in gene_symbol:
								gene_symbol.append(mytranscript.additional_annotations[0])


						#outline to use to output RUDs
						outline = str(entrez_id) + "\t" + chrom + "\t" + gene.strand + "\t" + str(basic_RUD) + "\t" + ",".join(map(str, RUDs)) + "\n"
						
						#outline to use to output polyA information for a species
						#outline = str(entrez_id) + "\t" + chrom + "\t" + gene.strand + "\t" + str(UTRbeginning) + "\t" + ",".join(map(str, PolyAsites)) + "\n"
					
						outf.write(outline)
	outf.close()		

def CUTR_vs_AUTR(true3UTRstarts, true3UTRends, reads, PorM, threshold):

	if plus.match(PorM):
		#get_RUDs_P requires true3UTRends to be sorted
		RUDs, basic_RUD, PolyAsites = get_RUDs_P(true3UTRends, reads, threshold)
		
	elif minus.match(PorM):
		#get_RUDs_N requires true3UTRstarts to be inversely sorted
		RUDs, basic_RUD, PolyAsites = get_RUDs_N(true3UTRstarts, reads, threshold)
		
	else:
		print "error in CUTR_vs_AUTR"

	return RUDs, basic_RUD, PolyAsites


def get_RUDs_P(ends, reads, threshold) :
	if len(ends) == 0:
		return [-9], -9, ends
	if len(ends) == 1:
		return [-7], -7, ends
	
	elif len(reads) < threshold:
		return [-1], -1, ends
	else:
		
		UTR_regions = []

		for i in ends:
			UTR_regions.append(0)

		for i in reads:
			k = 0
			for j in ends:
				if i <= j:
					UTR_regions[k] = UTR_regions[k] + 1
					break
				k = k + 1
		
		if  UTR_regions[0] == 0:
			return [-2], -2, ends
		else:
			cUTR = UTR_regions[0]
			
			del UTR_regions[0]
			
			RUD = []
			aUTRcount = 0
			for i in UTR_regions:
				RUD.append(i/float(cUTR))
				aUTRcount = aUTRcount + i
			basic_RUD = float(aUTRcount)/cUTR
			
			return RUD, basic_RUD, ends

#expects things to be reversed
def get_RUDs_N(starts, reads, threshold) :
	if len(starts) == 0:
		return [-9], -9, starts
	if len(starts) == 1:
		return [-7], -7, starts
	
	elif len(reads) < threshold:
		return [-1], -1, starts
	else:
		
		UTR_regions = []

		for i in starts:
			UTR_regions.append(0)

		for i in reads:
			k = 0
			for j in starts:
				if i >= j:
					UTR_regions[k] = UTR_regions[k] + 1
					break
				k = k + 1
		
		if  UTR_regions[0] == 0:
			return [-2], -2, starts
		else:
			
			cUTR = UTR_regions[0]
			
			del UTR_regions[0]
			
			RUD = []
			aUTRcount = 0
			for i in UTR_regions:
				RUD.append(i/float(cUTR))
				aUTRcount = aUTRcount + i
			basic_RUD = float(aUTRcount)/cUTR
			
			return RUD, basic_RUD, starts



def Mastertuplemaker(three_UTRs,listofvalues,chrom,PorM, extension):
	if plus.match(PorM):
		start, end = tuplechangerP(three_UTRs, extension)
		true3UTRstarts, true3UTRends, UTRregion_start, UTRregion_end = tupletesterP(start,end,listofvalues,chrom)
		UTRbeginning = UTRregion_start

	elif minus.match(PorM):
		start, end = tuplechangerN(three_UTRs, extension)
		true3UTRstarts, true3UTRends, UTRregion_start, UTRregion_end = tupletesterN(start,end,listofvalues,chrom)
		UTRbeginning = UTRregion_end

	else:
		print "Error in Mastertuplemaker"
		

	return true3UTRstarts, true3UTRends, UTRregion_start, UTRregion_end, UTRbeginning

def tuplechangerP(tuples, extension):
	tuple1 = []

	for i, j in tuples:
		tuple1.append(int(i))
		tuple1.append(j+extension)

	return tuple1[0],tuple1[1]

def tuplechangerN(tuples, extension):
	tuple1 = []

	for i, j in tuples:
		tuple1.append(i-extension)
		tuple1.append(int(j))

	return tuple1[0],tuple1[1]


def tupletesterP(start, end, listofvalues, chrom):
	list1 = []
	list2 = []
	
	for k, l in listofvalues:
		if k == chrom:
			z = int(l)				
			if z > start and z < end:
				list1.append(start)
				list2.append(z)
	
	if len(list2) > 0:
		list2.sort()
		UTRend = list2[-1]
	else:
		
		#probably want to change
		UTRend = end

	return list1, list2, start, UTRend

def tupletesterN(start, end, listofvalues, chrom):
	list1 = []
	list2 = []

	for k, l in listofvalues:
		if k == chrom:
			z = int(l)
			if z > start and z < end:
				list1.append(z)
				list2.append(end)
	
	if len(list1) > 0:				
		list1.sort(reverse = True)
		UTRstart = list1[-1]
	else:
		#probably want to change
		UTRstart = start

	return list1, list2, UTRstart, end


def associate_tags_with_3UTR (tag_positions, UTRregion_start, UTRregion_end):
	#Cannot use similar code from Utility_Extended as it requires strings while we are dealing with integers
	my_tag_list = []
	if (Utility_extended.is_list_sorted(tag_positions)==0):
		my_tag_list = sorted(tag_positions)
	else:
		my_tag_list = tag_positions
		
	assert (UTRregion_start<=UTRregion_end)
	start_ind = bisect.bisect_left(my_tag_list, UTRregion_start)
	end_ind = bisect.bisect_right(my_tag_list, UTRregion_end)
	tags = my_tag_list[start_ind : end_ind]
	
	return tags

def main(argv):
	parser = OptionParser()
	parser.add_option("-f", "--forwardreadfile", action="store", type="string", dest="ReadsOnForwardStrand", help="input bed file for RNASeq raw reads on forward strand", metavar="<file>")
	parser.add_option("-r", "--reversereadfile", action="store", type="string", dest="ReadsOnReverseStrand", help="input bed file for RNASeq raw reads on reverse strand", metavar="<file>")
	parser.add_option("-u", "--entrez_genes_file", action="store", type="string", dest="entrez_genes", metavar="<file>", help="file with curated known genes clustered by entrez ID in pickle format")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="outfile name", metavar="<file>")
	parser.add_option("-s", "--species", action="store", type="string", dest="species",help="species, mm8, hg18, etc", metavar="<str>")
	parser.add_option("-p", "--PAfile", action="store", type="string", dest="PAfile", help="input bed3 file", metavar="<file>")	
	parser.add_option("-e", "--extension", action="store", type="int", dest="extension",help="integer value denoting how far downstream the program should look for polyadenylation sites past the Entrez given 3'UTR end", metavar="<float>")
		

	(opt, args) = parser.parse_args(argv)

	if len(argv) < 14:
		parser.print_help()
		sys.exit(1)

	startTime = time.time()

	allowance = 10

	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species]
		chrom_lengths = GenomeData.species_chrom_lengths[opt.species]
	else:
		print "This species is not recognized, exiting"
		sys.exit(1)

	# entrez_gene_collection is a KnownEntrezGenes class object. The core is a entrez_genes.entrez_genes is a dic (keyed by entrez_id) of lists of EntrezGene object
	annotation = open(opt.entrez_genes, 'rb')
	entrez_gene_collection = Entrez.KnownEntrezGenes(chroms, pickle.load(annotation))
	annotation.close()

	# test module
	test = 0
	if test == 1:
		print "Testing gene structure"
		test_id = 79947
		Entrez.test_gene_structure(entrez_gene_collection, test_id)


	# Filter cluster of refseq_ids (keyed by entrez_id) according to the criterion of identical cdsEnd
	entrez_ids_with_unique_cdsEnd = entrez_gene_collection.get_ids_with_unique_cdsEnd()
	print "There are ", len(entrez_ids_with_unique_cdsEnd), " Entrez IDs each of which has a unique cdsEnd."


	#get total read count
	totalcount_F = get_total_tag_counts.get_total_tag_counts(opt.ReadsOnForwardStrand)
	totalcount_R = get_total_tag_counts.get_total_tag_counts(opt.ReadsOnReverseStrand)
	totalcount = totalcount_F + totalcount_R
	print totalcount_F, totalcount_R

	#Clear the file and write the first line
	outf = open(opt.outfile, 'w')
	
	#outline to use to output polyA information for a species	
	#outline = "# Entrez ID" + "\t" + "Chrom" + "\t" + "Strand" + "\t" + "UTRstart" + "\t" + "PolyAsites" + "\n"
	#outline to use to output RUDs
	outline = "# Entrez ID" + "\t" + "Chrom" + "\t" + "Strand" + "\t" + "Basic_RUD" + "\t" + "List_of_subRUDs" + "\n"
	outf.write(outline)
	outf.close()

	#index: column in bed file for sorting
	index = 2

	print "Process genes on forward strand"
	entrez_ids_on_forward_strand = entrez_gene_collection.get_strand_specific_ids("+", entrez_ids_with_unique_cdsEnd)
	print "There are ", len(entrez_ids_on_forward_strand), " Entrez IDs on forward strand."
	entrez_gene_subset = Entrez.KnownEntrezGenes(chroms, entrez_gene_collection.subset(entrez_ids_on_forward_strand))

	Calculate3UTRUsage(entrez_gene_subset, opt.ReadsOnForwardStrand, chroms, opt.outfile, allowance, opt.PAfile, opt.extension, index)


	print "Process genes on reverse strand"
	entrez_ids_on_reverse_strand = entrez_gene_collection.get_strand_specific_ids("-", entrez_ids_with_unique_cdsEnd)
	print "There are ", len(entrez_ids_on_reverse_strand), " Entrez IDs on reverse strand."
	entrez_gene_subset = Entrez.KnownEntrezGenes(chroms, entrez_gene_collection.subset(entrez_ids_on_reverse_strand))

	Calculate3UTRUsage(entrez_gene_subset, opt.ReadsOnReverseStrand, chroms, opt.outfile, allowance, opt.PAfile, opt.extension, index)

	print "it took", time.time() - startTime, "seconds."

if __name__ == "__main__":
	main(sys.argv)