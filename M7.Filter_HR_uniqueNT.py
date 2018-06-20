#!/usr/bin/python3

#quit if old python version
import sys
if sys.version_info[:2] < (3,0):
	print("This script is written for python 3.0 or above, please upgrade to a more recent version of python")
	exit()

#get command line arguments
import sys
INPUT = sys.argv[1]

#####
#set these variables to assign the core sequence and location to be used to search for HR/NOHR events
#sequences used are GFPStop-TTATAA; GFPSilent-CGCGCG; GFPWT-CGCGCC
coreStart=83 #59+24
coreEnd=88 #64+24
#HRmutations=("C40T","G41T","C42A","G43T","C44A","C45A")
refNTs=list("CGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCcgcgccGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGAC")
#####

#import re module for regular expressions and compile regex's
import re
Cregex=re.compile("\d+?\D{1}")

#set counters and process variables
totalseqs=MutTargetCounter=NOMUT=UNIQUEseqs=UNIQUEseqs2=UNIQUEoutwseqs=NOTUNIQUEseqs=NOTBUT=0
coreEndRange=coreEnd+1

#open OUTFILES
FILENAME1=INPUT+"unique.HR"
FILENAME2=INPUT+".moremuts.HR"
with open(FILENAME1,"w") as UNIQUE_OUT:
	with open(FILENAME2,"w") as NOTUNIQUE_OUT:
			#open SAMFILE
			with open(INPUT,"r") as SAMFILE:
				print(INPUT)
				for line in SAMFILE:
					if not line.startswith("M"):
						continue #skip SAM header lines
					totalseqs += 1
					SAMfields=line.split("\t") #parse SAM fields
					refstart=int(SAMfields[3])
					CIGAR=SAMfields[5]
					readNTs=list(SAMfields[9])
					readstart=1
					MutCounter = []
					MutCollector = []
					for j in range(0,len(refNTs)):
						if j >= coreStart-1:
							if j <= coreEnd:
								continue
						if readNTs[j] != refNTs[j]:
							MutCounter.append(j)
							MutCollector.append(readNTs[j])
					MutTargetCounter=UniqueMuts=0
					if (len(MutCounter)==0 or all (i == 'N' for i in MutCollector)):
						NOMUT += 1
						#print(MutCollector)
					elif len(MutCounter)==1: #
						for i in MutCounter:
							if 69 <= i <= 100:
								UNIQUEseqs += 1
								UNIQUE_OUT.write(line)
							else:
								UNIQUEoutwseqs += 1					
					else:
						for i in MutCounter:
							if 69 <= i <= 100:
								UniqueMuts += 1
						if UniqueMuts == 1:
							NOTBUT +=1
							UNIQUE_OUT.write(line)
						else: 
							NOTUNIQUEseqs += 1
							#print(MutCounter,MutCollector)
							NOTUNIQUE_OUT.write(line)
print(totalseqs,"mapped alignments processed\n",NOMUT,"mapped reads with no additional mutations.\n",UNIQUEseqs,"mapped reads with a unique mutation, written to",FILENAME1,"\n",UNIQUEoutwseqs,"mapped reads with a unique mutation outwith variable region. \n",NOTUNIQUEseqs,"mapped reads with multiple muts, written to ",FILENAME2,"\n",NOTBUT," mapped which have only a single mutation in the variable region.\n")
exit()

