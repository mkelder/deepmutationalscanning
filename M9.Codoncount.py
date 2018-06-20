#!/usr/bin/python3

#quit if old python version
import sys
if sys.version_info[:2] < (3,0):
	print("This script is written for python 3.0 or above, please upgrade to a more recent version of python")
	exit()

#get command line arguments
if len(sys.argv)==2:
	INPUT = sys.argv[1]
else:
	print("This python script counts nucleotide mutation frequencies, nucleotide mutations and amino acid mutations for all alignments in a SAM file relative to either the STOP or SILENT GFP sequences \n\nUsage:\npython mutation_SAM.py <SAMFILE> <STOP/SILENT>\n")
	exit()
if INPUT in ("-h", "--help"):
	print("This python script counts nucleotide mutation frequencies, nucleotide mutations and amino acid mutations for all alignments in a SAM file relative to either the STOP or SILENT GFP sequences \n\nUsage:\npython mutation_SAM.py <SAMFILE> <STOP/SILENT>\n")
#if "STOP" in INPUT2:
	#HR="Stop"
#elif "SILENT" in INPUT2:
	#HR="Silent"
#else:
	#print("This python script counts nucleotide mutation frequencies, nucleotide mutations and amino acid mutations for all alignments in a SAM file relative to either the STOP or SILENT GFP sequences \n\nUsage:\npython mutation_SAM.py <SAMFILE> <STOP/SILENT>\n")
	#exit()

#####
#set these variables to assign the core sequence and location to be used to search for HR/NOHR events
#sequences used are GFPStop-TTATAA; GFPSilent-CGCGCG; GFPWT-CGCGCC, trigger to use Stop or Silent sequences conferred from command line argument
WT_NTseq="CGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCcgcgccGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGAC"
#ReadingFrame is nucleotide frame (1,2,3) to translate WT_NTseq into amino acids
ReadingFrame=2
#outputstart is the nucleotide in WT_NTseq from which the output should start reporting (inclusive), set at 0 to inactivate and output from the start of WT_NTseq
outputstart=0
#outputstop is the nucleotide in WT_NTseq at which the output should stop reporting (inclusive), set at 0 to inactivate and output to the end of WT_NTseq
outputstop=0
#seqstart is the start nucleotide of the WT_NTseq sequence relative to a larger reference sequence (e.g. GFP ORF) for reporting nucleotide and amino acid mutations
seqstart=54
#ntstart is the codon at which it should start (for beta-catenin 13)
ntstart=13
####
#codons_twist ={'ATG':'M','AAC':'N','ATC':'I','CTG':'L','GCT':'A','TTC':'F','CAG':'Q','GGT':'G','GAT':'D','AAG':'K','ACC':'T','TGG':'W','AGA':'R','AGC':'S','TGC':'C','GTG':'V','CCT':'P','CAC':'H','GAG':'E','TAC':'Y'}


#define function to translate dna to protein
def translate_dna(dna_sequence,frame=None):
	if frame is None:
		frame=1
	codon_table ={'TCA':'S','TCC':'S','TCG':'S','TCT':'S','TTC':'F','TTT':'F','TTA':'L','TTG':'L','TAC':'Y','TAT':'Y','TAA':'*','TAG':'*','TGC':'C','TGT':'C','TGA':'*','TGG':'W','CTA':'L','CTC':'L','CTG':'L','CTT':'L','CCA':'P','CCC':'P','CCG':'P','CCT':'P','CAC':'H','CAT':'H','CAA':'Q','CAG':'Q','CGA':'R','CGC':'R','CGG':'R','CGT':'R','ATA':'I','ATC':'I','ATT':'I','ATG':'M','ACA':'T','ACC':'T','ACG':'T','ACT':'T','AAC':'N','AAT':'N','AAA':'K','AAG':'K','AGC':'S','AGT':'S','AGA':'R','AGG':'R','GTA':'V','GTC':'V','GTG':'V','GTT':'V','GCA':'A','GCC':'A','GCG':'A','GCT':'A','GAC':'D','GAT':'D','GAA':'E','GAG':'E','GGA':'G','GGC':'G','GGG':'G','GGT':'G'}
	protein_sequence=""
	for k in range(frame-1,len(dna_sequence),3):
		codon=dna_sequence[k:k+3]
		codon=codon.upper()
		if codon in codon_table:
			protein_sequence += codon_table[codon]
		else: 
			if "+" in codon:
				protein_sequence += ">"
				break
			elif "-" in codon:
				protein_sequence += "<"
				break  
			else:
				protein_sequence += "."
	return protein_sequence

#import re module for regular expressions and compile regex's
import re
Cregex=re.compile("\d+?\D{1}")

#set counters and process variables
#if HR=="Silent":
#	print("Using Silent codon")
#	HRNTmuts=set(["C333G"])
#	HRAAmuts=set(["A111A"])
#elif HR=="Stop":
#	print("Using Stop codon")
#	HRNTmuts=set(["C328T","G329T","C330A","G331T","C332A","C333A"])
#	HRAAmuts=set(["R110L","A111*"])
#else:
#	print("Error setting HR mutation type, please edit python script and follow instructions in line 16 to set HR mutation type")
#	quit() 
adjustNT=seqstart-1 #to adjust NT numbering for reporting mutations
codon_number=(seqstart+ReadingFrame+2-1) // 3
adjustAA=codon_number-1 #to adjust AA numbering for reporting mutations

if outputstart > len(WT_NTseq):
	outputstart = 1
elif outputstart == 0:
	outputstart = 1
outputstop=outputstop+1
if outputstop > len(WT_NTseq):
	outputstop = len(WT_NTseq)+1
elif outputstop <= outputstart:
	outputstop = len(WT_NTseq)+1
elif outputstop == 0:
	outputstop = len(WT_NTseq)+1

WT_AAseq=translate_dna(WT_NTseq,ReadingFrame)
WT_NTseq="."+WT_NTseq.upper() #. placed in position 0 of sequence to better align NT numbering with python index
WT_AAseq="."+WT_AAseq.upper() #. placed in position 0 of sequence to better align AA numbering with python index
endSeqRange=len(WT_NTseq)
WTAAs={}
WTcodons={}
TCA={}
TCC={}
TCG={}
TCT={}
TTC={}
TTT={}
TTA={}
TTG={}
TAC={}
TAT={}
TAA={}
TAG={}
TGC={}
TGT={}
TGA={}
TGG={}
CTA={}
CTC={}
CTG={}
CTT={}
CCA={}
CCC={}
CCG={}
CCT={}
CAC={}
CAT={}
CAA={}
CAG={}
CGA={}
CGC={}
CGG={}
CGT={}
ATA={}
ATC={}
ATT={}
ATG={}
ACA={}
ACC={}
ACG={}
ACT={}
AAC={}
AAT={}
AAA={}
AAG={}
AGC={}
AGT={}
AGA={}
AGG={}
GTA={}
GTC={}
GTG={}
GTT={}
GCA={}
GCC={}
GCG={}
GCT={}
GAC={}
GAT={}
GAA={}
GAG={}
GGA={}
GGC={}
GGG={}
GGT={}
AAcons={}
AAdel={}
AAnons={}
AAother={}
total={}

#A={}
#C={}
#G={}
#T={}
#d={}
#i={}
#fs={}
#sil={}
#mis={}
#non={}

for j in range (ntstart,33):
	WTAAs[j]=int(0)
	WTcodons[j]=int(0)
	TCA[j]=int(0)
	TCC[j]=int(0)
	TCG[j]=int(0)
	TCT[j]=int(0)
	TTC[j]=int(0)
	TTT[j]=int(0)
	TTA[j]=int(0)
	TTG[j]=int(0)
	TAC[j]=int(0)
	TAT[j]=int(0)
	TAA[j]=int(0)
	TAG[j]=int(0)
	TGC[j]=int(0)
	TGT[j]=int(0)
	TGA[j]=int(0)
	TGG[j]=int(0)
	CTA[j]=int(0)
	CTC[j]=int(0)
	CTG[j]=int(0)
	CTT[j]=int(0)
	CCA[j]=int(0)
	CCC[j]=int(0)
	CCG[j]=int(0)
	CCT[j]=int(0)
	CAC[j]=int(0)
	CAT[j]=int(0)
	CAA[j]=int(0)
	CAG[j]=int(0)
	CGA[j]=int(0)
	CGC[j]=int(0)
	CGG[j]=int(0)
	CGT[j]=int(0)
	ATA[j]=int(0)
	ATC[j]=int(0)
	ATT[j]=int(0)
	ATG[j]=int(0)
	ACA[j]=int(0)
	ACC[j]=int(0)
	ACG[j]=int(0)
	ACT[j]=int(0)
	AAC[j]=int(0)
	AAT[j]=int(0)
	AAA[j]=int(0)
	AAG[j]=int(0)
	AGC[j]=int(0)
	AGT[j]=int(0)
	AGA[j]=int(0)
	AGG[j]=int(0)
	GTA[j]=int(0)
	GTC[j]=int(0)
	GTG[j]=int(0)
	GTT[j]=int(0)
	GCA[j]=int(0)
	GCC[j]=int(0)
	GCG[j]=int(0)
	GCT[j]=int(0)
	GAC[j]=int(0)
	GAT[j]=int(0)
	GAA[j]=int(0)
	GAG[j]=int(0)
	GGA[j]=int(0)
	GGC[j]=int(0)
	GGG[j]=int(0)
	GGT[j]=int(0)
	AAcons[j]=int(0)#accounted for in line 248
	AAdel[j]=int(0)
	AAnons[j]=int(0)
	AAother[j]=int(0)	
	total[j]=int(0)

for j in range(ntstart,33): #Setup reference codons
	total[j] += 1
	k=3*j-2+ReadingFrame-1
	wtNT=WT_NTseq[k]+WT_NTseq[k+1]+WT_NTseq[k+2]
	WT_AA=translate_dna(wtNT,1)
	WTcodons[j]=wtNT
	WTAAs[j]=WT_AA

#open OUTFILES
FILENAME1=INPUT+".codoncount"
with open(FILENAME1,"w") as COUNT:
	COUNT.write("\t\t\t\tS\tS\tS\tS\tF\tF\tL\tL\tY\tY\t*\t*\tC\tC\t*\tW\tL\tL\tL\tL\tP\tP\tP\tP\tH\tH\tQ\tQ\tR\tR\tR\tR\tI\tI\tI\tM\tT\tT\tT\tT\tN\tN\tK\tK\tS\tS\tR\tR\tV\tV\tV\tV\tA\tA\tA\tA\tD\tD\tE\tE\tG\tG\tG\tG\n")
	COUNT.write("Pos\tWT_AA\tWT_codon\tTotal_count\tTCA\tTCC\tTCG\tTCT\tTTC\tTTT\tTTA\tTTG\tTAC\tTAT\tTAA\tTAG\tTGC\tTGT\tTGA\tTGG\tCTA\tCTC\tCTG\tCTT\tCCA\tCCC\tCCG\tCCT\tCAC\tCAT\tCAA\tCAG\tCGA\tCGC\tCGG\tCGT\tATA\tATC\tATT\tATG\tACA\tACC\tACG\tACT\tAAC\tAAT\tAAA\tAAG\tAGC\tAGT\tAGA\tAGG\tGTA\tGTC\tGTG\tGTT\tGCA\tGCC\tGCG\tGCT\tGAC\tGAT\tGAA\tGAG\tGGA\tGGC\tGGG\tGGT\n")
	with open(INPUT,"r") as SAMFILE:
		print("Reading",INPUT)
		for line in SAMFILE:
			if line.startswith("@"):
				continue #skip SAM header lines
			mutations=[]
			AAmutations=[]
			SAMfields=line.split("\t") #parse SAM fields
			ID=SAMfields[0]
			refstart=int(SAMfields[3])
			CIGAR=SAMfields[5]
			readNTs=list(SAMfields[9])
			readstart=1
			seq={}
			for j in range(0,endSeqRange): #initialise dictionary for nt positions, 0 used to bin any trailing read nucleotides
				seq[j]="."
			CIGARlist=Cregex.findall(CIGAR) #parse CIGAR string
			for CIGARpart in CIGARlist:
				bases=str()
				for char in CIGARpart:
					if char.isdigit():
						bases += str(char)
					else :
						code=str(char)
				refstart=1 #############################################
				refend=refstart+int(bases)
				readend=readstart+int(bases)
				if code == "D": #deal with deletions in read
					for j in range(refstart,refend):
						seq[j]="-"
						#d[j] += 1 
						mut=WT_NTseq[j]+str(j+adjustNT)+"-"
						mutations.append(mut)
				elif code == "S": #discard soft clipped nucleotides in read
					refend=refstart
					for j in range(readstart,readend):
						readNTs.pop(0)
				elif code == "I": #for insertions, mark preceding nucleotide with +, and concatenate with lower case inserted nucleotides N+ins
					refend=refstart
					nt=refstart-1
					seq[nt]+="+"
					#i[nt] += 1
					mut=WT_NTseq[j]+str(j+adjustNT)+"+"
					for j in range(readstart,readend):
						insN=readNTs.pop(0)
						mut += insN.lower()
					mutations.append(mut)
				elif code == "M": #skip through 'matched' nucleotides in read
					for j in range(refstart,refend):
						seq[j]=readNTs.pop(0)
				else :
					print("\nError, unrecognised code in CIGAR string:",code,"\n")
					exit()
				refstart=refend
				readstart=readend
			read_NTseq=""
			for j in range(ntstart,33): #now that insertions/deletions are dealt with, compare through each triplet(element in list) with its aligned WT nucleotide (character in string)	
				total[j] += 1
				k=3*j-2+ReadingFrame-1
				mNT=seq[k]+seq[k+1]+seq[k+2]
				wtNT=WT_NTseq[k]+WT_NTseq[k+1]+WT_NTseq[k+2]				
				if mNT==wtNT: #no change
					AAcons[j] += 1
				elif "TCA"  in mNT: 
					TCA[j] += 1
				elif "TCC"  in mNT: 
					TCC[j] += 1
				elif "TCG"  in mNT: 
					TCG[j] += 1
				elif "TCT"  in mNT: 
					TCT[j] += 1
				elif "TTC"  in mNT: 
					TTC[j] += 1
				elif "TTT"  in mNT: 
					TTT[j] += 1
				elif "TTA"  in mNT: 
					TTA[j] += 1
				elif "TTG"  in mNT: 
					TTG[j] += 1
				elif "TAC"  in mNT: 
					TAC[j] += 1
				elif "TAT"  in mNT: 
					TAT[j] += 1
				elif "TAA"  in mNT: 
					TAA[j] += 1
				elif "TAG"  in mNT: 
					TAG[j] += 1
				elif "TGC"  in mNT: 
					TGC[j] += 1
				elif "TGT"  in mNT: 
					TGT[j] += 1
				elif "TGA"  in mNT: 
					TGA[j] += 1
				elif "TGG"  in mNT: 
					TGG[j] += 1
				elif "CTA"  in mNT: 
					CTA[j] += 1
				elif "CTC"  in mNT: 
					CTC[j] += 1
				elif "CTG"  in mNT: 
					CTG[j] += 1
				elif "CTT"  in mNT: 
					CTT[j] += 1
				elif "CCA"  in mNT: 
					CCA[j] += 1
				elif "CCC"  in mNT: 
					CCC[j] += 1
				elif "CCG"  in mNT: 
					CCG[j] += 1
				elif "CCT"  in mNT: 
					CCT[j] += 1
				elif "CAC"  in mNT: 
					CAC[j] += 1
				elif "CAT"  in mNT: 
					CAT[j] += 1
				elif "CAA"  in mNT: 
					CAA[j] += 1
				elif "CAG"  in mNT: 
					CAG[j] += 1
				elif "CGA"  in mNT: 
					CGA[j] += 1
				elif "CGC"  in mNT: 
					CGC[j] += 1
				elif "CGG"  in mNT: 
					CGG[j] += 1
				elif "CGT"  in mNT: 
					CGT[j] += 1
				elif "ATA"  in mNT: 
					ATA[j] += 1
				elif "ATC"  in mNT: 
					ATC[j] += 1
				elif "ATT"  in mNT: 
					ATT[j] += 1
				elif "ATG"  in mNT: 
					ATG[j] += 1
				elif "ACA"  in mNT: 
					ACA[j] += 1
				elif "ACC"  in mNT: 
					ACC[j] += 1
				elif "ACG"  in mNT: 
					ACG[j] += 1
				elif "ACT"  in mNT: 
					ACT[j] += 1
				elif "AAC"  in mNT: 
					AAC[j] += 1
				elif "AAT"  in mNT: 
					AAT[j] += 1
				elif "AAA"  in mNT: 
					AAA[j] += 1
				elif "AAG"  in mNT: 
					AAG[j] += 1
				elif "AGC"  in mNT: 
					AGC[j] += 1
				elif "AGT"  in mNT: 
					AGT[j] += 1
				elif "AGA"  in mNT: 
					AGA[j] += 1
				elif "AGG"  in mNT: 
					AGG[j] += 1
				elif "GTA"  in mNT: 
					GTA[j] += 1
				elif "GTC"  in mNT: 
					GTC[j] += 1
				elif "GTG"  in mNT: 
					GTG[j] += 1
				elif "GTT"  in mNT: 
					GTT[j] += 1
				elif "GCA"  in mNT: 
					GCA[j] += 1
				elif "GCC"  in mNT: 
					GCC[j] += 1
				elif "GCG"  in mNT: 
					GCG[j] += 1
				elif "GCT"  in mNT: 
					GCT[j] += 1
				elif "GAC"  in mNT: 
					GAC[j] += 1
				elif "GAT"  in mNT: 
					GAT[j] += 1
				elif "GAA"  in mNT: 
					GAA[j] += 1
				elif "GAG"  in mNT: 
					GAG[j] += 1
				elif "GGA"  in mNT: 
					GGA[j] += 1
				elif "GGC"  in mNT: 
					GGC[j] += 1
				elif "GGG"  in mNT: 
					GGG[j] += 1
				elif "GGT"  in mNT: 
					GGT[j] += 1
				elif "ATG"  in mNT: 
					ATG[j] += 1
				elif "AAC"  in mNT: 
					AAC[j] += 1
				elif "ATC"  in mNT: 
					ATC[j] += 1
				elif "CTG"  in mNT: 
					CTG[j] += 1
				elif "GCT"  in mNT: 
					GCT[j] += 1
				elif "TTC"  in mNT: 
					TTC[j] += 1
				elif "CAG"  in mNT: 
					CAG[j] += 1
				elif "GGT"  in mNT: 
					GGT[j] += 1
				elif "GAT"  in mNT: 
					GAT[j] += 1
				elif "AAG"  in mNT: 
					AAG[j] += 1
				elif "ACC"  in mNT: 
					ACC[j] += 1
				elif "TGG"  in mNT: 
					TGG[j] += 1
				elif "AGA"  in mNT: 
					AGA[j] += 1
				elif "AGC"  in mNT: 
					AGC[j] += 1
				elif "TGC"  in mNT: 
					TGC[j] += 1
				elif "GTG"  in mNT: 
					GTG[j] += 1
				elif "CCT"  in mNT: 
					CCT[j] += 1
				elif "CAC"  in mNT: 
					CAC[j] += 1
				elif "GAG"  in mNT: 
					GAG[j] += 1
				elif "TAC"  in mNT: 
					TAC[j] += 1
				elif "TAC"  in mNT: 
					TAC[j] += 1		
				elif "." in mNT:		
					AAdel[j] += 1	
				elif "TAA" or "TAG" or "TGA" in mNT:		
					AAnons[j] += 1	
				else:		
					AAother[j] += 1	
					print(AAother)
		for j in range(ntstart,33):
			#total[j]=A[j]+C[j]+G[j]+T[j]+d[j]+i[j]
			if total[j]==0:
				total[j]=1
			AAcount=j+82
			perc_TCA=TCA[j]/float(total[j])
			perc_TCC=TCC[j]/float(total[j])
			perc_TCG=TCG[j]/float(total[j])
			perc_TCT=TCT[j]/float(total[j])
			perc_TTC=TTC[j]/float(total[j])
			perc_TTT=TTT[j]/float(total[j])
			perc_TTA=TTA[j]/float(total[j])
			perc_TTG=TTG[j]/float(total[j])
			perc_TAC=TAC[j]/float(total[j])
			perc_TAT=TAT[j]/float(total[j])
			perc_TAA=TAA[j]/float(total[j])
			perc_TAG=TAG[j]/float(total[j])
			perc_TGC=TGC[j]/float(total[j])
			perc_TGT=TGT[j]/float(total[j])
			perc_TGA=TGA[j]/float(total[j])
			perc_TGG=TGG[j]/float(total[j])
			perc_CTA=CTA[j]/float(total[j])
			perc_CTC=CTC[j]/float(total[j])
			perc_CTG=CTG[j]/float(total[j])
			perc_CTT=CTT[j]/float(total[j])
			perc_CCA=CCA[j]/float(total[j])
			perc_CCC=CCC[j]/float(total[j])
			perc_CCG=CCG[j]/float(total[j])
			perc_CCT=CCT[j]/float(total[j])
			perc_CAC=CAC[j]/float(total[j])
			perc_CAT=CAT[j]/float(total[j])
			perc_CAA=CAA[j]/float(total[j])
			perc_CAG=CAG[j]/float(total[j])
			perc_CGA=CGA[j]/float(total[j])
			perc_CGC=CGC[j]/float(total[j])
			perc_CGG=CGG[j]/float(total[j])
			perc_CGT=CGT[j]/float(total[j])
			perc_ATA=ATA[j]/float(total[j])
			perc_ATC=ATC[j]/float(total[j])
			perc_ATT=ATT[j]/float(total[j])
			perc_ATG=ATG[j]/float(total[j])
			perc_ACA=ACA[j]/float(total[j])
			perc_ACC=ACC[j]/float(total[j])
			perc_ACG=ACG[j]/float(total[j])
			perc_ACT=ACT[j]/float(total[j])
			perc_AAC=AAC[j]/float(total[j])
			perc_AAT=AAT[j]/float(total[j])
			perc_AAA=AAA[j]/float(total[j])
			perc_AAG=AAG[j]/float(total[j])
			perc_AGC=AGC[j]/float(total[j])
			perc_AGT=AGT[j]/float(total[j])
			perc_AGA=AGA[j]/float(total[j])
			perc_AGG=AGG[j]/float(total[j])
			perc_GTA=GTA[j]/float(total[j])
			perc_GTC=GTC[j]/float(total[j])
			perc_GTG=GTG[j]/float(total[j])
			perc_GTT=GTT[j]/float(total[j])
			perc_GCA=GCA[j]/float(total[j])
			perc_GCC=GCC[j]/float(total[j])
			perc_GCG=GCG[j]/float(total[j])
			perc_GCT=GCT[j]/float(total[j])
			perc_GAC=GAC[j]/float(total[j])
			perc_GAT=GAT[j]/float(total[j])
			perc_GAA=GAA[j]/float(total[j])
			perc_GAG=GAG[j]/float(total[j])
			perc_GGA=GGA[j]/float(total[j])
			perc_GGC=GGC[j]/float(total[j])
			perc_GGG=GGG[j]/float(total[j])
			perc_GGT=GGT[j]/float(total[j])
			perc_conc=AAcons[j]/float(total[j])
			perc_AAdel=AAdel[j]/float(total[j])
			perc_AAnons=AAnons[j]/float(total[j])
			perc_AAother=AAother[j]/float(total[j])			
			COUNT.write("%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\n" %(AAcount,WTAAs[j],WTcodons[j],total[j],perc_TCA,perc_TCC,perc_TCG,perc_TCT,perc_TTC,perc_TTT,perc_TTA,perc_TTG,perc_TAC,perc_TAT,perc_TAA,perc_TAG,perc_TGC,perc_TGT,perc_TGA,perc_TGG,perc_CTA,perc_CTC,perc_CTG,perc_CTT,perc_CCA,perc_CCC,perc_CCG,perc_CCT,perc_CAC,perc_CAT,perc_CAA,perc_CAG,perc_CGA,perc_CGC,perc_CGG,perc_CGT,perc_ATA,perc_ATC,perc_ATT,perc_ATG,perc_ACA,perc_ACC,perc_ACG,perc_ACT,perc_AAC,perc_AAT,perc_AAA,perc_AAG,perc_AGC,perc_AGT,perc_AGA,perc_AGG,perc_GTA,perc_GTC,perc_GTG,perc_GTT,perc_GCA,perc_GCC,perc_GCG,perc_GCT,perc_GAC,perc_GAT,perc_GAA,perc_GAG,perc_GGA,perc_GGC,perc_GGG,perc_GGT,perc_conc,perc_AAdel,perc_AAnons,perc_AAother))
print("\nDone.Amino acid mutation frequencies written to",FILENAME1,"\n\n")
exit()
