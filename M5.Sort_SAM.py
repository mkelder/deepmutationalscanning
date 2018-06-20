#!/usr/bin/python3

#quit if old python version
import sys
if sys.version_info[:2] < (3,0):
    print("This script is written for python 3.0 or above, please upgrade to a more recent version of python")
    exit()

#get command line arguments
import sys
if len(sys.argv)==4:
    INPUT = sys.argv[1]
    INPUT2 = sys.argv[2]
else:
    print("This python script sorts SAM files depending on whether the alignments posses a wild-type, homologous repair, or ambiguous 6nt core sequence. Specify whether to use STOP or SILENT HR sequence in the second command line argument\n\nUsage:\npython Sort_SAM_HR.py <SAMFILE> <STOP/SILENT>\n")
    exit()
if INPUT in ("-h", "--help"):
    print("This python script sorts SAM files depending on whether the alignments posses a wild-type, homologous repair, or ambiguous 6nt core sequence. Specify whether to use STOP or SILENT HR sequence in the second command line argument\n\nUsage:\npython Sort_SAM_HR.py <SAMFILE> <STOP/SILENT>\n")
if "STOP" in INPUT2:
    HR="Stop"
elif "SILENT" in INPUT2:
    HR="Silent"
else:
    print("This python script counts nucleotide mutation frequencies, nucleotide mutations and amino acid mutations for all alignments in a SAM file relative to either the STOP or SILENT GFP sequences \n\nUsage:\npython mutation_SAM.py <SAMFILE> <STOP/SILENT>\n")
    exit()

#####
#set these variables to assign the core sequence and location to be used to search for HR/NOHR events
#sequences used are GFPStop-TTATAA; GFPSilent-CGCGCG; GFPWT-CGCGCC
if "STOP" in INPUT2:
    HRcore="TTATAA"
elif "SILENT" in INPUT2:
    HRcore="CGCGCG"
NOHRcore="CGCGCC"
coreStart=83 #GFP default: 83
coreEnd=88 #GFP default: 88
WTsequence=sys.argv[3]
#####

#import re module for regular expressions and compile regex's
import re
Cregex=re.compile("\d+?\D{1}")

#set counters and process variables
totalseqs=HRseqs=NOHRseqs=AMBseqs=0
coreEndRange=coreEnd+1

#open OUTFILES
FILENAME1=INPUT+"."+INPUT2+".HR"
FILENAME2=INPUT+"."+INPUT2+".NOHR"
FILENAME3=INPUT+"."+INPUT2+".AMB"
with open(FILENAME1,"w") as HR_OUT:
    with open(FILENAME2,"w") as NOHR_OUT:
        with open(FILENAME3,"w") as AMB_OUT:
            #open SAMFILE
            with open(INPUT,"r") as SAMFILE:
                print("Reading",INPUT)
                for line in SAMFILE:
                    if not line.startswith("M"):
                        continue #skip SAM header lines
                    totalseqs += 1
                    SAMfields=line.split("\t") #parse SAM fields
                    refstart=int(SAMfields[3])
                    CIGAR=SAMfields[5]
                    readNTs=list(SAMfields[9])
                    readstart=1
                    seq={}
                    for j in range(coreStart,coreEndRange):
                        seq[j]="."
                    CIGARlist=Cregex.findall(CIGAR)
                    for CIGARpart in CIGARlist:
                        bases=str()
                        for char in CIGARpart:
                            if char.isdigit():
                                bases += str(char)
                            else :
                                code=str(char)
                        refend=refstart+int(bases)
                        readend=readstart+int(bases)
                        if code == "D":
                            for j in range(refstart,refend):
                                if coreStart<=j<=coreEnd:
                                    seq[j]="-" 
                        elif code == "S":
                            refend=refstart
                            for j in range(readstart,readend):
                                readNTs.pop(0)
                        elif code == "I":
                            refend=refstart
                            nt=refstart-1
                            if coreStart<=nt<=coreEnd:
                                seq[nt]+="+"
                            for j in range(readstart,readend):
                                readNTs.pop(0)
                        elif code == "M":
                            for j in range(refstart,refend):
                                mNT=readNTs.pop(0)
                                if coreStart<=j<=coreEnd:
                                    seq[j]=mNT
                        else :
                            print("\nError, unrecognised code in CIGAR string:",code,"\n")
                            exit()
                        refstart=refend
                        readstart=readend
                    core=""
                    for j in range(coreStart,coreEndRange):
                        core+=str(seq[j])
                    if core==HRcore:
                        HRseqs += 1
                        HR_OUT.write(line)
                    elif core==NOHRcore:
                        NOHRseqs += 1
                        NOHR_OUT.write(line)
                    elif core=="------":
                        NOHRseqs += 1
                        NOHR_OUT.write(line)
                    else :
                        AMBseqs += 1
                        AMB_OUT.write(line)
print("Done\n\n",totalseqs,"alignments processed\n",HRseqs,"mapped reads with HR core, written to",FILENAME1,"\n",NOHRseqs,"mapped reads with NOHR (WT) core, written to",FILENAME2,"\n",AMBseqs,"mapped reads with ambiguous core, written to",FILENAME3)
exit()

