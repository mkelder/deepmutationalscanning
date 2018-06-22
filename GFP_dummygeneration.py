#!/usr/bin/python3
#This python script generates a series of mock FASTQ sequencing reads from a reference FASTA sequence. 
#The original sequence will be outputted 10 times then for each base, 10 sequences more sequences will be outputted with a random nucleotide substitution at that base. 
#Each sequences will have randomised per nucleotide quality scores >= 20 and be outputted as paired end sequences to dummy_seqs_1.fastq and dummy_seqs_2.fastq"


#import modules
import sys
import random
import re

#quit if old python version
if sys.version_info[:2] < (3,0):
	print("This script is written for python 3.0 or above, please upgrade to a more recent version of python")
	exit()
	
INFOTEXT="\nThis python script generates a series of mock FASTQ sequencing reads from a reference FASTA sequence. The original sequence will be outputted 10 times then for each base, 10 sequences more sequences will be outputted with a random nucleotide substitution at that base. Each sequences will have randomised per nucleotide quality scores >= 20 and be outputted as paired end sequences to dummy_seqs_1.fastq and dummy_seqs_2.fastq\n"
USAGE="\nUsage:\npython dummy_mutant_sequences.py INPUTSEQUENCE.fasta\n"

#set label variable to use for sequence names
LABEL="WTGFP"
#set length of sequencing reads required
READLENGTH=150
#randomly select quality scores from this set	
QUALS=["e","f","g","h","i","j"]
#represent each fusion sequence this number of times in the output file
REPS=1
#Set how many mutations we want to introduce
MUTNUMBER=15
#define valid regex for FASTA inputs
VALIDFASTA=re.compile("[^ACGTacgt]")
#set types of mutations ("Sub","Ins","Del")
MUTCLASS=["Sub","Ins","Del"]	
#set number of reads per mutation number
READNUMBER=200

#define reverse complementation function
def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]

#get command line arguments for input file (FASTA)
if len(sys.argv)==2:
	INPUTFILE1=sys.argv[1]
else:
	print("%s%s" %(INFOTEXT,USAGE))
	exit()
if INPUTFILE1 in ("-h", "--help"):
	print("%s%s" %(INFOTEXT,USAGE))
	exit()

#check INPUT file exists
try:
	open(INPUTFILE1,"r")
except IOError:
	print("\nERROR: Invalid input file: ",INPUTFILE1, " cannot be accessed\n")
	print("%s" %(USAGE))
	exit()

#check that temporary OUTPUT file can be written
for m in range(0,MUTNUMBER+1):
	WRITETEMPFILE="temporary_working_sequences"+str(m)+str(".fasta")
	try:
		open(WRITETEMPFILE,"w")
	except IOError:
		print("\nERROR: Cannot write to temporary file: ",WRITETEMPFILE, "\n")
		print("%s" %(USAGE))
		exit()

#check that final OUTPUT file can be written
OUTPUT1="dummy_all_muts_readerrors_1.fastq"
OUTPUT2="dummy_all_muts_readerrors_2.fastq"
try:
	open(OUTPUT1,"w")
except IOError:
	print("\nERROR: Invalid output file: ",OUTPUT1, " cannot be accessed\n")
	print("%s" %(USAGE))
	exit()
try:
	open(OUTPUT2,"w")
except IOError:
	print("\nERROR: Invalid output file: ",OUTPUT2, " cannot be accessed\n")
	print("%s" %(USAGE))
	exit()

#load INPUTFILE1, check FASTA format, and process into a single string
TESTFIRSTLINE=1
INSEQ1=str("")
with open(INPUTFILE1,"r") as READFILE:
	for line in READFILE:
		if TESTFIRSTLINE==1:
			if not line.startswith(">"):
				print("\nERROR: Invalid input file: ",INPUTFILE1, " not FASTA format\n")
				print("%s" %(USAGE))
				exit()
			TESTFIRSTLINE=0
		else:
			line=line.rstrip()
			if VALIDFASTA.search(line):
				print("\nERROR: Invalid character: ", VALIDFASTA.findall(line), "in FASTA file: ",INPUTFILE1, "\n")
				print("\nERROR: Invalid input file: ",INPUTFILE1, " not FASTA format\n")
				print("%s" %(USAGE))
				exit()
			else:
				INSEQ1 += str(line)
READFILE.close()

#determine length of reference sequence
REFLENGTH=len(INSEQ1)
if len(INSEQ1) == 0:
	print("\nERROR: input file: ",INPUTFILE1, " has no sequence information\n")
	print("%s%s" %(INFOTEXT,USAGE))
	exit()
if READLENGTH > REFLENGTH:
	READLENGTH = REFLENGTH
#convert reference sequence to all lowercase
REFSEQ=INSEQ1.lower()

#Create 20000 dummy PCR sequences with MUTNUMBER of substitutions and output these to temporary fasta files
for z in MUTCLASS:
	for m in range(0,MUTNUMBER+1):
		if m==0:
			WRITETEMPFILE="temporary_working_sequences_"+str(z)+str("_")+str(m)+str(".fasta")
			with open (WRITETEMPFILE,"w") as TEMP2:
				MUTPOS=0
				SEQNAME1=str(LABEL)+str("_")+str(z)+str("_")+str(m)+str("_")+str(MUTPOS)
				TEMP2.write(">%s\n%s\n" %(SEQNAME1,INSEQ1))
			TEMP2.close()
		else:
			WRITETEMPFILE="temporary_working_sequences"+str("_")+str(z)+str(m)+str(".fasta")
			with open (WRITETEMPFILE,"w") as TEMP:
				if z == "Sub":
					for i in range(0,READNUMBER+1):
						MUTSEQ1=INSEQ1.lower()
						target=random.sample(range(1,157),m)
						target.sort()
						MUTPOS=(str(target)[1:-1]).replace(" ","")
						for n in target:
							TARGETBASE=MUTSEQ1[n]
							MUTANTBASE=random.choice(['A','C','G','T'])
							while (MUTANTBASE.lower() == TARGETBASE):
									MUTANTBASE=random.choice(['A','C','G','T'])
							MUTSEQ1=str(MUTSEQ1[:n])+str(MUTANTBASE)+str(INSEQ1[n+1:])
						SEQNAME1=str(LABEL)+str("_")+str(z)+str("_")+str(m)+str("_")+MUTPOS+str("_")+str(i)
						TEMP.write(">%s\n%s\n" %(SEQNAME1,MUTSEQ1))
				elif z == "Del":
					for i in range(0,READNUMBER+1):
						MUTSEQ1=INSEQ1.lower()
						target=random.sample(range(1,157),m)
						target.sort()
						MUTPOS=(str(target)[1:-1]).replace(" ","")
						c = 0
						for n in target:
							TARGETBASE=MUTSEQ1[n-c]
							c += 1
							MUTANTBASE=random.choice(['A','C','G','T'])
							MUTSEQ1=str(MUTSEQ1[:n+1])+str(MUTANTBASE)+str(MUTSEQ1[n+1:])
						SEQNAME1=str(LABEL)+str("_")+str(z)+str("_")+str(m)+str("_")+MUTPOS+str("_")+str(i)
						TEMP.write(">%s\n%s\n" %(SEQNAME1,MUTSEQ1))
				elif z == "Ins":
					for i in range(0,READNUMBER+1):
						MUTSEQ1=INSEQ1.lower()
						target=random.sample(range(1,157),m)
						target.sort()
						MUTPOS=(str(target)[1:-1]).replace(" ","")
						c = 0
						for n in target:
							TARGETBASE=MUTSEQ1[n+c]
							c += 1
							MUTANTBASE=random.choice(['A','C','G','T'])
							MUTSEQ1=str(MUTSEQ1[:n+1])+str(MUTANTBASE)+str(MUTSEQ1[n+1:])
						SEQNAME1=str(LABEL)+str("_")+str(z)+str("_")+str(m)+str("_")+MUTPOS+str("_")+str(i)
						TEMP.write(">%s\n%s\n" %(SEQNAME1,MUTSEQ1))
			TEMP.close()
	
#Use the generated PCR dummy sequences as temples and generate REPS number of dummy sequence reads, and output these to fastq files.
with open(OUTPUT1,"w") as OUTFILE1:
	with open(OUTPUT2,"w") as OUTFILE2:
		for x in MUTCLASS:
			for m in range(1,MUTNUMBER+1):
				READTEMPFILE="temporary_working_sequences_"+str(x)+str(m)+str(".fasta")
				with open (READTEMPFILE,"r") as TEMPFILE:
					for line in TEMPFILE:
						if line.startswith(">"):
							line=line.rstrip()
							SEQNAME=line
						else:
							line=line.rstrip()
							SEQ=line
							for i in range (0,REPS):
								SEQNAME1=SEQNAME+str("_Rep")+str(i)
								SEQ1=SEQ[:READLENGTH]
								QUALSCORES1=""
								for k in range (0,READLENGTH):
									QUALSCORES1 += str(random.choice(QUALS))
								SEQNAME2=SEQNAME1
								SEQ2=revcomp(SEQ[-READLENGTH:])
								QUALSCORES2=""
								for k in range (0,READLENGTH):
									QUALSCORES2 += str(random.choice(QUALS))
								OUTFILE1.write("@%s\n%s\n+%s\n%s\n" %(SEQNAME1,SEQ1,SEQNAME1,QUALSCORES1))
								OUTFILE2.write("@%s\n%s\n+%s\n%s\n" %(SEQNAME2,SEQ2,SEQNAME2,QUALSCORES2))
					TEMPFILE.close()
		OUTFILE1.close()
		OUTFILE2.close()
exit()
