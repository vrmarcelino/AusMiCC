#! /usr/bin/python3

from Bio import SeqIO
import argparse
import re
import subprocess
import os
import shutil
import ntpath
import re
from pathlib import Path
from ete3 import NCBITaxa

# local imports
from ausmicc_scripts import fconnector
from ausmicc_scripts import path_structure
from ausmicc_scripts import ctables
from ausmicc_scripts import fparsedb
from ausmicc_scripts import fsample
from ausmicc_scripts import fisolate
from ausmicc_scripts import fadd16S

import time
start_time = time.time()

from Bio import SeqIO

__version__ = "version 2.0"

databases=path_structure.ref_db()
blastdb = databases.blastdb_16Smicro

ncount = 1
cutoff = ['-g']
sequenceHash = {}
sequenceHash_single_end = {}
outputfile = "AllSequences.fna"


#Cutoff Hash Values
cutoffhash = {}
cutoffhash[1] = ["-l", "0.2", "-u", "0.5"]
cutoffhash[2] = ["-l", "0.5", "-u", "1.0"]
cutoffhash[3] = ["-l", "1.0", "-u", "1.5"]
cutoffhash[4] = ["-l", "1.3", "-u", "2"]
cutoffhash[5] = ["-l", "1.5", "-u", "2.5"]
cutoffhash[6] = ["-l", "2.2", "-u", "3.5"]
cutoffhash[7] = ["-l", "2.6", "-u", "4"]
cutoffhash[8] = ["-l", "3", "-u", "4"]
cutoffhash[9] = ["-l", "3.4", "-u", "4.5"]
cutoffhash[10] = ["-l", "4", "-u", "5"]
cutoffhash[11] = ["-l", "4.3", "-u", "5.5"]
cutoffhash[12] = ["-l", "4.7", "-u", "6"]
cutoffhash[13] = ["-l", "5.1", "-u", "6.5"]
cutoffhash[14] = ["-l", "5.5", "-u", "7"]
cutoffhash[15] = ["-l", "6.0", "-u", "7.5"]
cutoffhash[16] = ["-l", "6.4", "-u", "8"]
cutoffhash[17] = ["-l", "6.8", "-u", "8.5"]
cutoffhash[18] = ["-l", "7.2", "-u", "9"]
cutoffhash[19] = ["-l", "7.6", "-u", "9.5"]
cutoffhash[20] = ["-l", "8", "-u", "10"]

#These values are used to keep track of which sequences have been added with only 7f or 1510 to ensure sequences are combined where possible
FwdIdHash = {}
RevIdHash = {}
CombinedIdHash = {}

oldnames = {}
newnames = {}
replacementnames = {}


parser = argparse.ArgumentParser()
parser.add_argument('-d', dest='directory', help='A directory to process and import all sequences into single file')
parser.add_argument('-f', dest='file', help='A file to import sequences from')
parser.add_argument('-t', dest='threads', type=int, default=4, help='The number threads to use for analysis [4]')
parser.add_argument('-q', dest='quality', type=int, default=40, help='Average sequence quality score below which a sequence will be filtered. [40]')
parser.add_argument('-p', dest='paired',  action='store_true', help='Only keep that have both forward and reverse pairs')
parser.add_argument('-m', dest='minlength', default=400, type=int, help='The minimum length sequence to retain [400]')
parser.add_argument('-w', dest='windowsize', default=5, type=int, help='The window size for N removal [5]')
parser.add_argument('-M', dest='match', action='store_true', help='Enforce matching of known 16S rRNA region')
parser.add_argument('-n', dest='ncount', type=int, default=1, help='The number of Ns in a given window size [-w] where trimming is performed [1]')
parser.add_argument('-b', dest='blastdb', help='The path to the blast database to use for analysis [' + blastdb + ']')
parser.add_argument('-c', dest='cutoff', type=int, help='Percentage cutoff for OTU assignment (Must be 1, 2, 3, 5, 8 or 10) [5]')
parser.add_argument('-k', dest='known', help='A file containing a single column of the sequence ids where the isolate has already been stored and had DNA extracted. This list will be used to populate the known OTU columns in the Results.txt file')
parser.add_argument('-B', dest='bootstrap', type=int, default=100, help='Number of bootstraps to use when generating phylogenetic tree. [100]')
parser.add_argument('-H', dest='leaveheader', action='store_true', help='Use this flag to specify not to alter the headers other than removing _7f/_1510r')
parser.add_argument('-S', dest='modifyfile', help='Tab delimited file in the format: <CurrentName>\t<NewName>')
parser.add_argument('-F', dest='forward', action='store_true', help='Use this flag to add _7f to each sequence automatically')
parser.add_argument('-R', dest='reverse', action='store_true', help='Use this flag to add _1510r to each sequence automatically')

parser.add_argument('--taxonomy', default='Y', help="Use this flag to include the taxonomy step in analysis, options are 'Y' and 'N'. default=Y", required=False)
parser.add_argument('--add2db', default='Y', help="Add data to the AusMiCC database, options are 'Y' and 'N'. default=Y", required=False)

parser.add_argument('--align', dest='align', action='store_true', help='Use this flag to include alignment step in analysis')
parser.add_argument('--tree', dest='tree', action='store_true', help='Use this flag to include tree generation step in analysis')

parser.add_argument('-v', dest='verbose', action='store_true', default = False, help='Print more messages of progress and status. Useful for debugging.')

parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__))

args = parser.parse_args()

if args.directory is None and args.file is None:
	print("Neither directory (-d) or file (-f) containing fasta sequences is defined. One must be defined")
	exit(-1)
	
if args.directory is None:
	outputfile = args.file

if args.cutoff is None:
	args.cutoff = 5
	
if args.cutoff is not None and args.cutoff not in cutoffhash:
	print("Acceptable cutoff values are between 1 and 20 (integers only)")
	print("Please check values are re-run")
	exit(-1)

if args.forward and args.reverse:
	print("Please specify either -F or -R not both")
	exit(-1)
	
outputfile = ntpath.basename(outputfile)
outputfile = os.path.splitext(outputfile)[0]
outputfile = outputfile + ".fna"


#Functions
def modify_file_clean_name(name):
	name = name.replace("\n", "")
	name = name.replace(" ", "")
	return name	

def modify_file_import():
	mf = open(args.modifyfile, "r")
	for line in mf.readlines():
		line = line.replace("\n", "")
		a = line.split("\t")
		on = modify_file_clean_name(a[0])
		nn = modify_file_clean_name(a[1])
		if args.verbose is not None: print("Old Name: " + a[0] + " New Name: " + a[1])
	
		if on in oldnames:
			print("ERROR: Replacement name contains duplicates. Please check " + args.modifyfile + " and rerun analysis")
			print("\tLine: " + line)
			exit(-1)
		oldnames[on] = 0
		
		if nn in newnames:
			print("ERROR: Replacement name contains duplicates. Please check " + args.modifyfile + " and rerun analysis")
			print("\tLine: " + line)
			exit(-1)
		newnames[nn] = 0
		replacementnames[on] = nn

def modify_header(header):
	if replacementnames is not None:
		header_seqname = header
		if os.sep in header_seqname:
			header_seqname = header_seqname[header_seqname.index(os.sep)+1:]
		if args.verbose: print("Searching Header: " + header_seqname)
		if header_seqname in replacementnames:
			if args.verbose: print("Replacing Header: " + header)
			header = header.replace(header_seqname, replacementnames[header_seqname])
			if args.verbose: print("Replacing with: " + header)
		else:
			print("WARNING: Header not found in replacement file: " + header)
	
	else:
		if args.forward is not None:
			if "_7f" not in header[-3:].lower():
				header = header + "_7f"
		if args.reverse is not None:
			if "_1510r" not in header[-6:].lower():
				header = header + "_1510r"
	return header

def check_header(header):
	if not isForward(header) and not isReverse(header):
		print("WARNING: Header must specify _7f or _1510r. Skipping: " + header)
		return False
	return True		

def clean_header(header):
	header = header.replace(">", "")
	header = header.replace(".ab1", "")
	header = header.replace(".txt", "")
	header = header.replace("\n", "")
	header = header.replace("_7f", "")
	header = header.replace("_7F", "")
	header = header.replace("_1510r", "")
	header = header.replace("_1510R", "")
	if not args.leaveheader:
		if "-" in header:
			header = header[header.index("-"):]
		if "_" in header:
			header = header[header.index("_")+1:]
		if "_" in header:
			header = header[header.index("_")+1:]
		a = header.split('\t')
		if len(a) > 1:
			header = a[0]	
	return header

#function for removing N's from the start of the sequence (start = True or start = False)
def remove_Ns(seq, start):
	if not start:
		#reverse the sequence
		seq = seq[::-1]

	trimming = True
	i = 0
	invalidcharacters = ['N', 'n', 'W', 'w', 'M', 'm']

	while args.windowsize < len(seq) and trimming:
		count = 0
		for character in invalidcharacters:
			count = count + seq.count(character,0, args.windowsize)
		if count >= args.ncount:
			seq = seq[1:]
		else:
			trimming = False
			
	if not start:
		#reverse the sequence
		seq = seq[::-1]

	return seq
		

def clean_sequence(header, seq, forward):
	#First check it only contains valid characters
	if not valid_seq_characters(seq):
		if args.verbose:
			print(header + ": Invalid characters")
		return None
		
	#First remove N's
	seq = remove_Ns(seq, True)
	seq = remove_Ns(seq, False)

	#Check the length is ok after removing sequences
	if len(seq) < args.minlength:
		if args.verbose:
			print("Removing short sequence:  " + header + " " + str(len(seq)))
		return None
	
	return seq		

#Function to create temporary fasta file based on sequence
def create_fasta_file(filepath, header, sequence):
	f1 = open(filepath, "w")
	f1.write(">" + header + "\n")
	f1.write(sequence)
	f1.close()

#Function to check if header is forward sequence
def isForward(header):
	if "_7F" in header or "_7f" in header:
		return True
	else:
		return False

#Function to check if header is reverse sequence
def isReverse(header):
	if "_1510R" in header or "_1510r" in header:
		return True
	else:
		return False

#Function to check sequence only contains characters defined in regex search (ATCGatgcNn)
def valid_seq_characters(strg, search=re.compile(r'[^ATCGatcgNnMmWw]').search):
	return not bool(search(strg))

#Function responsible for adding sequences to sequenceHash. This function calls the relevant QC methods
def add_sequence(header, sequence):
	if check_header(header):
		if args.modifyfile is not None:
			header = modify_header(header)
		if args.verbose is not None: print("\tAdding Sequence: " + header)
		sequence = clean_sequence(header, sequence, isForward(header))
		if sequence is not None:
			if args.verbose is not None: print("\tAdding Sequence of length: " + str(len(sequence))) 
			id = clean_header(header)
			#If the ID has not been seen before
			if id not in sequenceHash:
				if isForward(header):
					FwdIdHash[id] = 1
					header_note = "_Fwd"
				if isReverse(header):
					RevIdHash[id] = 1
					header_note = "_Rev"
				sequenceHash[id] = sequence
				# also store in a dict that will keep the single-end sequences
				id_for_single_end = id + header_note
				sequenceHash_single_end[id_for_single_end] = sequence
			else:
				if id not in CombinedIdHash:
					#if the value being added is the same orientation as those that are already present just check length
					if (isForward(header) and id in FwdIdHash) or (isReverse(header) and id in RevIdHash):
						if len(sequenceHash[id]) < len(sequence):
							sequenceHash[id] = sequence
						else:
							if args.verbose is not None:
								print("Skipping " + id + " as shorter than matching sequencing") 
					#If this is different to the current header
					else:
						#If this is the reverse sequence create the f.fna file from the existing sequence
						if isReverse(header):
							create_fasta_file("f.fna", id, sequenceHash[id])
							create_fasta_file("r.fna", id, sequence)
							header_note = "_Rev"
						else:
							create_fasta_file("f.fna", id, sequence)
							create_fasta_file("r.fna", id, sequenceHash[id])
							header_note = "_Fwd"

						print("Merging: " + header)
						p = subprocess.Popen(["merger", "f.fna", "r.fna", "-sreverse2", "-outseq", "merged.fna", "-outfile", "temp.merger", "-gapopen", "100", "-gapextend", "2"], stdout=subprocess.PIPE)
						p.wait()
						#Replace sequence record with merged file
						for seq_record in SeqIO.parse("merged.fna", "fasta"):
							sequenceHash[id] = str(seq_record.seq.upper())
				
						#Store the fact sequences were combined
						CombinedIdHash[id] = 1

						# add info to single-end dictionary:
						id_for_single_end = id + header_note
						sequenceHash_single_end[id_for_single_end] = sequence

						#Remove temporary files
						os.remove("f.fna")
						os.remove("r.fna")
						os.remove("temp.merger")
						os.remove("merged.fna")
				else:
					if args.verbose:
						print("Skipping " + id + "(Not Sure why)")

	
#Check file type
def is_file_type(filepath, extensions):
	filename = ntpath.basename(filepath)
	if os.path.splitext(filename)[1][1:] in extensions:
		return True
	else:
		return False

#Check if this is a fasta file type
def is_fasta_sequence_filetype(filename):
	extensions = ["fna", "fa", "fasta", "txt"]
	return is_file_type(filename, extensions)
	
#Check if this is a seq type. The .seq format does not have a fasta header.
def is_seq_filetype(filename):
	extensions = ["seq"]
	return is_file_type(filename, extensions)

#Check if this is a ab1 type. As of version 1-15 this is the only sequence format accepted.
def is_ab1_sequence_filetype(filename):
	extensions = ["ab1"]
	return is_file_type(filename, extensions)

#Function to populate sequence hash from file of sequences
def import_sequence_from_ab1(file):
	record=SeqIO.read(file,"abi-trim")
	q = record.format("qual")

	total = 0
	bases = 0

	for line in q.split("\n"):
		if ">" not in line:
			line = line.replace("\n", "")
			val = line.split(" ")
			for v in val:
				if v != "":
					total = total + int(v)
					bases = bases + 1

	if total/bases > args.quality:
		add_sequence(file.replace(".ab1", ""), str(record.seq))
		return file
	else:
		if args.verbose:
			print("Skipping " + record.id + " due to low sequence quality: " + str(total/bases))


#Generic method for importing sequence file
def import_sequence_file(filename):
	if is_ab1_sequence_filetype(filename):
		if args.verbose : print("Importing ab1 File")
		added_files = import_sequence_from_ab1(filename)
		return added_files
	if is_fasta_sequence_filetype(filename):
		print("WARNING: Skipping " + filename + " as it is not an ab1 file")
	if is_seq_filetype(filename):
		print("WARNING: Skipping " + filename + " as it is not an ab1 file")


		
#Generic method for safe file delete
def delete_file(filename):
	if os.path.exists(filename):
		os.remove(filename)

# Function that converts a strain taxid to a species-level taxid.
# also returns the species-level classification:
def get_sp_taxid(query_taxid):
	if query_taxid != 'NOVEL':
		ncbi = NCBITaxa()
		query_taxid = int(query_taxid)
		lineage = ncbi.get_lineage(query_taxid)
		ranks = ncbi.get_rank(lineage)
		names = ncbi.get_taxid_translator(lineage)

		for key, val in ranks.items():
			if val == 'species':
				sp_taxid = key
				sp_name = names[key]
	else:
		sp_taxid = 0
		sp_name = 'unk_sp'

	species_info = [sp_taxid,sp_name]
	return species_info



########################################################################################################
#     Import the files, ensure they are reasonable quality and write cleaned fasta file
########################################################################################################
print("Importing Files....")
if args.verbose : print("\t--- %s seconds ---" % (time.time() - start_time))
section_start_time = time.time()


#If Modify File has been specified set this up
if args.modifyfile is not None:
	if args.verbose: print("Modifying file name headers using " + args.modifyfile)
	modify_file_import()

#Import sequences and ensure they are reasonable quality
imported_files = []
if args.file is not None:
	try_import = import_sequence_file(args.file)
	if try_import != None:
		imported_files.append(try_import)

if args.directory is not None:
	for filename in os.listdir(args.directory):
		if args.verbose is not None: print(filename)
		try_import = import_sequence_file(args.directory + os.sep + filename)
		if try_import != None:
			imported_files.append(try_import)

#if paired has been specified filter file
if args.paired:
	keys = list()
	for k in sequenceHash:
		keys.append(k)
	for k in keys:
		if k not in CombinedIdHash:
			sequenceHash.pop(k, None)
			if args.verbose:
				missingpair = "Reverse"
				if k in RevIdHash:
					missingpair = "Forward"	
				print("Removing " + k + " due to missing " + missingpair + " pair")
	

#Write the cleaned fasta file
cleanedfasta = outputfile
cleanedfasta = cleanedfasta.replace(".fasta", ".fna")
cleanedfasta = cleanedfasta.replace(".fa", ".fna")
cleanedfasta = cleanedfasta.replace(".fna", ".cleaned.fna")


cleanedfastaoutput = open(cleanedfasta, "w")
for key in sequenceHash:
	print (key)
	cleanedfastaoutput.write(">" + key + "\n" + sequenceHash[key] + "\n")
cleanedfastaoutput.close()

print("Finished importing sequences.")
if args.verbose : print("\t--- %s Section seconds ---" % (time.time() - section_start_time))

print("\t" + str(len(sequenceHash)) + " sequences imported")

if len(sequenceHash) == 0:
	print("Error: No sequences imported")
	exit(-1)
	
#################################################
#     Compare Sequences and determine OTUs
#################################################

print("Calculating OTUs...")
section_start_time = time.time()
if args.verbose : print("\t--- %s seconds ---" % (time.time() - start_time))


#The cutoff is -g unless args.cutoff has been specified


cropargs = ["CROP", "-i", cleanedfasta]

if len(sequenceHash) > 75:
	cropargs = cropargs + ["-b", str(int(len(sequenceHash)/50))]
	
if args.cutoff is not None and args.cutoff in cutoffhash:
	cropargs = cropargs + cutoffhash[args.cutoff]
		
if args.verbose is not None:
	print("Arguments: " + str(cropargs))

p = subprocess.Popen(cropargs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
p.wait()

delete_file("LikelihoodRatio.txt")
delete_file(cleanedfasta + ".unique")
delete_file(cleanedfasta + ".unique.list")
delete_file(cleanedfasta + ".unique.TempCenters.Rare")
if args.verbose : print("\t--- %s Section seconds ---" % (time.time() - section_start_time))


######################################
#     Blast Sequence
######################################
print("Blasting Sequences...", flush=True)
if args.verbose : print("\t--- %s seconds ---" % (time.time() - start_time))
section_start_time = time.time()

matchedspecies = open("MatchedSpecies.txt", "w")
matchedspecies.write("Genome ID\tBlast Match\tSequence Length\tAlignment Length\tIdentity\tE-value\tCoverage\n")
p = subprocess.Popen(["blastn", "-num_threads", str(args.threads), "-db", blastdb, "-query", cleanedfasta, "-outfmt", "6 qseqid stitle qlen length pident evalue qcovs", "-max_hsps", "1", "-max_target_seqs", "1"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
matchedspecies.write(p.stdout.read().decode("utf-8"))
matchedspecies.close()

#Create Results.txt File
#Populate blashhash with the blast results and add results to database
matchedspeciesin = open("MatchedSpecies.txt", "r")
blasthash = {}
for blastline in matchedspeciesin.readlines():
	if "Genome ID" not in blastline:
		blastline = blastline.replace("\n", "")
		a = blastline.split("\t")
		blasthash[a[0]] = blastline	
matchedspeciesin.close()

#If a cultured list has been provided populate the cultured list from this file
knownhash = {}
if args.known is not None:
	culturedin = open(args.known, "r")
	for knownline in culturedin.readlines():
		knownline = knownline.replace("\n", "")
		knownhash[knownline] = 1
		print("Adding '" + knownline + "' to knownhash", flush=True)


#Now write results
results = open("Results.txt", "w")
results.write("Genome ID\tBlast Match\tSequence Length\tAlignment Length\tIdentity\tE-value\tCoverage\tOTU Count\tCultured OTUs\tAll OTUs\n")


otufile = open(cleanedfasta + ".cluster.list", "r")


otucount = 0

for otuline in otufile.readlines():
	a = otuline.split("\t")
	
	if a[0] in blasthash:
		#Get the count of OTUs in this group
		a[1] = a[1].replace("\n", "")
		allotus = a[1].split(",")
		count = len(allotus)
		
		knownotus = ""
		#If provided determine if any members of this OTU have been cultured
		for o in allotus:
			if o in knownhash:
				if knownotus == "":
					knowotus = o
				else:
					knownotus = "," + knownotus
			
		
		results.write(blasthash[a[0]] + "\t" + str(count) + "\t" + knownotus + "\t" + a[1] + "\n")
		otucount = otucount + 1
	else:
		print("\t\tExcluding " + a[0] + " due to lack of blast hit", flush=True)
		
	
results.close()
otufile.close()

print("\tFound " + str(otucount) + " OTUs for analysis", flush=True)
if args.verbose : print("\t--- %s Section seconds ---" % (time.time() - section_start_time))


#if len(sequenceHash) == 1:
#	print("\tOnly 1 sequences provided. Skipping Alignment and Tree Building")
#	exit(-1)

######################################
#     Generate Taxonomy File
######################################
if args.taxonomy == 'Y':
	print("Generating Taxonomy File...", flush=True)
	if args.verbose : print("\t--- %s seconds ---" % (time.time() - start_time))
	section_start_time = time.time()
	p = subprocess.Popen(["taxonomy_gettaxon", "-r", cleanedfasta, "-o", cleanedfasta.replace(".fna", ".taxonomy")])
	p.wait()
	if args.verbose : print("\t--- %s Section seconds ---" % (time.time() - section_start_time))

elif args.taxonomy == 'N':
	pass
else:
	print ("--taxonomy must be either 'N' or 'Y'. Skipping taxids calculcation")
	
######################################
#     Perform alignment
######################################
if (args.align or args.tree) and len(sequenceHash) > 1:
	print("Performing alignment...", flush=True)
	if args.verbose : print("\t--- %s seconds ---" % (time.time() - start_time))
	section_start_time = time.time()
	#First remove previous ssuoutput if they exist
	if os.path.exists("ssuoutput"):
		print("\tRemoving previous ssuoutput")
		shutil.rmtree("ssuoutput")

	p = subprocess.Popen(["ssu-align", "-f", cleanedfasta, "ssuoutput"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	p.wait()

	alignment = cleanedfasta.replace(".fna", ".aln")
	alignmentoutput = open(alignment, "w")

	ssufile = open("ssuoutput/ssuoutput.bacteria.stk", "r")
	for line in ssufile.readlines():
		line = line.replace("\n", "")
		if not "#" in line[0:1] and line != "" and line != "//":
			a = line.split()
			a[1] = a[1].replace(".", "-")
			alignmentoutput.write(">" + a[0] + "\n" + a[1] + "\n")


	alignmentoutput.close()

	#First remove previous ssuoutput if they exist
	if os.path.exists("ssuoutput" + os.sep + "ssuoutput.nomatch"):
		print("\tUnable to find 16S rRNA gene in the following sequences:")
		nfh = open("ssuoutput" + os.sep + "ssuoutput.nomatch", "r")
		for nf in nfh.readlines():
			nf = nf.replace("\n", "")
			print("\t\t" + nf)

	if os.path.exists("ssuoutput") and args.verbose == False:
		shutil.rmtree("ssuoutput")
	
	if args.verbose : print("\t--- %s Section seconds ---" % (time.time() - section_start_time))


######################################
#     Create Tree
######################################
if args.tree:
	if len(sequenceHash) > 5:
		print("Creating tree...", flush=True)
		if args.verbose : print("\t--- %s seconds ---" % (time.time() - start_time))
		section_start_time = time.time()
		treename = alignment.replace(".aln", ".tree")
		if args.verbose:
			print("Running: " + "raxmlHPC-PTHREADS", "­­print­identical­sequences", "-p", "12345", "-T", str(args.threads), "-#", str(args.bootstrap), "-m", "GTRGAMMA", "-b", "12345", "-k", "-s", alignment, "-n", treename)
		p = subprocess.Popen(["raxmlHPC-PTHREADS", "­­print­identical­sequences", "-p", "12345", "-T", str(args.threads), "-#", str(args.bootstrap), "-m", "GTRGAMMA", "-b", "12345", "-k", "-s", alignment, "-n", treename], stdout=subprocess.PIPE)
		p.wait()

		if not args.verbose:
			os.rename("RAxML_bootstrap." + treename, treename)
			for p in Path(".").glob("RAxML_*"):
				p.unlink()
		if args.verbose : print("\t--- %s Section seconds ---" % (time.time() - section_start_time))

	else:
		print("Skipping tree generation due to lack of different species")

print("Completed Analysis")
if args.verbose : print("\t--- %s Seconds ---" % (time.time() - start_time))


######################################
#    Add info to AusMiCC database
######################################

if args.taxonomy == 'Y':
	print ("\nAdding info to the AusMiCC database\n")

	### connect to the database:
	aus_db_conn = fconnector.db_connection()
	cursor = aus_db_conn.cursor()

	# read taxonomy file and save taxids as a dictionary:
	taxonomy_fp = cleanedfasta.replace(".fna", ".taxonomy")

	names2tax_dic = {}
	with open(taxonomy_fp, 'r') as a:
		for line in a:
			if len(line) > 1:
				iso_name = line.split("\t")[0]
				iso_taxid = line.split("\t")[1]

				names2tax_dic[iso_name] = iso_taxid

	new_db_entries = []
	for in_file in imported_files:

		entry_inf = ctables.info_16s()

		if isForward(in_file):
			primer = "7f"
		else:
			primer = "1510r"

		#define db location, the sample name, copy file and store info in the entry_inf class obj:
		in_file_no_dir = in_file.split("/")[-1] # input file name without directory

		isolate_name_with_path = clean_header(in_file)
		isolate_name = isolate_name_with_path.split("/")[-1]

		amplicon_paths = path_structure.amplicon_paths()
		db_loc_ab1 = amplicon_paths.ab1_collec + "/" + in_file_no_dir

		entry_inf.isolate_name = isolate_name
		entry_inf.ab1_loc = db_loc_ab1
		entry_inf.primer = primer

		# rsync files:
		print ("rsync'ing files")
		p = subprocess.Popen(["rsync", "-avzh", in_file, db_loc_ab1], stdout=subprocess.PIPE)
		p.wait()
		print("saved %s file" % (in_file))

		# get partial and full len sequences:
		if primer == "7f":
			isolate_name_with_path_primer = isolate_name_with_path + "_Fwd"
		else:
			isolate_name_with_path_primer = isolate_name_with_path + "_Rev"
		entry_inf.sequence = sequenceHash_single_end[isolate_name_with_path_primer]

		#  check if it has the full len sequence:
		if isolate_name_with_path in CombinedIdHash:
			entry_inf.full_len_sequence = CombinedIdHash[isolate_name_with_path]

		# get taxids and species:
		lca_taxid = names2tax_dic[isolate_name_with_path]

		species_n_taxid = get_sp_taxid(lca_taxid)

		entry_inf.lca_taxid = lca_taxid
		entry_inf.sp_taxid = species_n_taxid[0]
		entry_inf.species = species_n_taxid[1]

		new_db_entries.append(entry_inf)

	# now add to the MySQL:
	print ("\nInteracting with MySQL\n")
	for obj in new_db_entries:
		print(obj.isolate_name, obj.species)

		# add info to the 16S table:
		fadd16S.add_16S_record(obj,cursor)

		#update the isolate table with species classification and full length sequences:
		fadd16S.update_isolate_info(obj,cursor)

		print ("added %s %s to the AusMiCC database" %(obj.isolate_name, obj.primer))

	aus_db_conn.commit()
	aus_db_conn.close()
	print ("\nDone!!\n")


#todo: if taxonomy == Y, maybe add a check in the begining to see if the isolate_names can be found in the db

