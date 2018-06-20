#!/usr/bin/python3

#quit if old python version
import sys
if sys.version_info[:2] < (3,0):
	print("This script is written for python 3.0 or above, please upgrade to a more recent version of python")
	exit()

#get command line arguments
import sys
if len(sys.argv)==2:
	INPUT = sys.argv[1]
else:
	print("This python script filters HR reads with only a single nt change. Specify whether to use STOP or SILENT HR sequence in the second command line argument\n\nUsage:\npython Sort_SAM_HR.py <SAMFILE> <STOP/SILENT>\n")
	exit()
if INPUT in ("-h", "--help"):
	print("This python script filters HR reads with only a single nt change. Specify whether to use STOP or SILENT HR sequence in the second command line argument\n\nUsage:\npython Sort_SAM_HR.py <SAMFILE> <STOP/SILENT>\n")

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
totalseqs=MutTargetCounter=NOMUTC=MUTC=0
coreEndRange=coreEnd+1

#open OUTFILES
FILENAME1=INPUT+"."+".withmutations.NOHR"
FILENAME2=INPUT+"."+".WT.NOHR"
with open(FILENAME1,"w") as MUT_OUT:
	with open(FILENAME2,"w") as WT_OUT:
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
					if len(MutCounter)==0:
						NOMUTC += 1
						WT_OUT.write(line)
						#print(MutCollector)				
					else:
						MUTC +=1
						MUT_OUT.write(line)
print(totalseqs,"mapped alignments processed\n",NOMUTC,"mapped reads with no additional mutations written to",FILENAME2,"\n")
exit()

