#!/usr/bin/python3

#quit if old python version
import sys
if sys.version_info[:2] < (3,0):
	print("This script is written for python 3.0 or above, please upgrade to a more recent version of python")
	exit()

#get command line arguments
if len(sys.argv)==3:
	INPUT = sys.argv[1]
	INPUT2 = sys.argv[2]
else:
	print("This python script counts nucleotide mutation frequencies, nucleotide mutations and amino acid mutations for all alignments in a SAM file relative to either the STOP or SILENT GFP sequences \n\nUsage:\npython mutation_SAM.py <SAMFILE> <STOP/SILENT>\n")
	exit()
if INPUT in ("-h", "--help"):
	print("This python script counts nucleotide mutation frequencies, nucleotide mutations and amino acid mutations for all alignments in a SAM file relative to either the STOP or SILENT GFP sequences \n\nUsage:\npython mutation_SAM.py <SAMFILE> <STOP/SILENT>\n")
if "STOP" in INPUT2:
	HR="Stop"
elif "SILENT" in INPUT2:
	HR="Silent"
else:
	print("This python script counts nucleotide mutation frequencies, nucleotide mutations and amino acid mutations for all alignments in a SAM file relative to either the STOP or SILENT GFP sequences \n\nUsage:\npython mutation_SAM.py <SAMFILE> <STOP/SILENT>\n")
	exit()

#####
#set these variables to assign the core sequence and location to be used to search for HR/NOHR events
#sequences used are GFPStop-TTATAA; GFPSilent-CGCGCG; GFPWT-CGCGCC, trigger to use Stop or Silent sequences conferred from command line argument
WT_NTseq="CGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCcgcgccGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGAC"
#STOP_NTseq="CGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCttataaGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGAC"
#ReadingFrame is nucleotide frame (1,2,3) to translate WT_NTseq into amino acids
ReadingFrame=2
#outputstart is the nucleotide in WT_NTseq from which the output should start reporting (inclusive), set at 0 to inactivate and output from the start of WT_NTseq
outputstart=36
#outputstop is the nucleotide in WT_NTseq at which the output should stop reporting (inclusive), set at 0 to inactivate and output to the end of WT_NTseq
outputstop=133
#seqstart is the start nucleotide of the WT_NTseq sequence relative to a larger reference sequence (e.g. GFP ORF) for reporting nucleotide and amino acid mutations
seqstart=246
####

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
if HR=="Silent":
	print("Using Silent codon")
	HRNTmuts=set(["C333G"])
	HRAAmuts=set(["A111A"])
elif HR=="Stop":
	print("Using Stop codon")
	HRNTmuts=set(["C328T","G329T","C330A","G331T","C332A","C333A"])
	HRAAmuts=set(["R110L","A111*"])
else:
	print("Error setting HR mutation type, please edit python script and follow instructions in line 16 to set HR mutation type")
	quit() 
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

nc={}
A={}
C={}
G={}
T={}
d={}
i={}
fs={}
sil={}
mis={}
non={}

for j in range (0,endSeqRange):
	nc[j]=int(0)
	A[j]=int(0)
	C[j]=int(0)
	G[j]=int(0)
	T[j]=int(0)
	d[j]=int(0)
	i[j]=int(0)
	fs[j]=int(0)
	sil[j]=int(0)
	mis[j]=int(0)
	non[j]=int(0)

#open OUTFILES
FILENAME1=INPUT+".NTcount"
FILENAME2=INPUT+".mutations"
with open(FILENAME1,"w") as COUNT:
	COUNT.write("pos\tWT_NT\tWT_AA\tA\tC\tG\tT\tdel\tins\tnc\ttotal\t%A\t%C\t%G\t%T\t%del\t%ins\t%nc\t%SILENT\t%MISSENSE\t%NONSENSE\t%INDEL\n")
	with open(FILENAME2,"w") as MUT:
		MUT.write("ID\tNT_MUT_WT\tAA_MUT_WT\tNT_MUT_HR\tAA_MUT_HR\n")
		#open SAMFILE
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
					refend=refstart+int(bases)
					readend=readstart+int(bases)
					if code == "D": #deal with deletions in read
						for j in range(refstart,refend):
							seq[j]="-"
							d[j] += 1 
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
						i[nt] += 1
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
				for j in range(1,endSeqRange): #now that insertions/deletions are dealt with, compare through each read nucleotide (element in list) with its aligned WT nucleotide (character in string)
					mNT=seq[j]
					if mNT==WT_NTseq[j]: #no change
						nc[j] += 1
					else : #identify and deal with substitutions
						if "."  not in mNT:
							if "+"  not in mNT:
								if "-"  not in mNT:
									mut=WT_NTseq[j]+str(j+adjustNT)+mNT
									mutations.append(mut)
					if mNT == "A": #count nucleotide frequency
						A[j] += 1
					elif mNT == "C":
						C[j] += 1
					elif mNT == "G":
						G[j] += 1
					elif mNT == "T":
						T[j] += 1
					if "+" in mNT: #remove +ins from read nucleotides where insertions have occured to allow downstream numbering to align with WT later during mutation analysis
						mNT=mNT.split("+")[0]
					read_NTseq += mNT #generate aligned read sequence
				read_NTseq = "."+read_NTseq #insert . at position 0 like for WT sequence
				for mut in mutations: #loop through list of NT mutations
					mut=mut[1:]
					NTpos=""
					NTcode=""
					for char in mut:
						if char.isdigit():
							NTpos += str(char)
						else :
							NTcode += str(char)
					NTpos=int(NTpos)-adjustNT
					if "+" in NTcode:
						NTpos += 1
					codon_number=1+((NTpos-ReadingFrame) // 3)
					codonStart=ReadingFrame+(3*(codon_number-1))
					nextCodonStart=codonStart+3
					readCodon=read_NTseq[codonStart:nextCodonStart] #get and translate read and WT codons containing mutated nucleotides
					WTCodon=WT_NTseq[codonStart:nextCodonStart]
					readAA=translate_dna(readCodon,1)
					WTAA=translate_dna(WTCodon,1)
					if "+" in NTcode: #deal with insertions (+ins removed from read NTs, record kept in mutation list)
						readAA = ">"
						NTpos -= 1
						fs[NTpos] += 1
					elif "-" in NTcode:
						readAA = "<"
						fs[NTpos] += 1
					else :
						if "." not in readAA:
							if ">" not in readAA:
								if "<" not in readAA:
									if "*" in readAA:
										non[NTpos] += 1
									elif WTAA == readAA:
										sil[NTpos] += 1
									else :
										mis[NTpos] += 1
					AAmut=WTAA+str(codon_number+adjustAA)+readAA
					AAmutations.append(AAmut) #generate AA mutations
				muts=",".join(set(mutations)) #write AA and NT mutation lists to file
				AAmuts=",".join(set(AAmutations))
				mutsHR=",".join(set(mutations).difference(HRNTmuts))
				AAmutsHR=",".join(set(AAmutations).difference(HRAAmuts))
				MUT.write("%s\t%s\t%s\t%s\t%s\n" % (ID,muts,AAmuts,mutsHR,AAmutsHR))
			total={} #prepare and write counts/NT mutation frequencies to file
			aa={}
			for j in range(outputstart,outputstop):
				total[j]=A[j]+C[j]+G[j]+T[j]+d[j]+i[j]
				if total[j]==0:
					total[j]=1
				percA=100*A[j]/float(total[j])
				percC=100*C[j]/float(total[j])
				percG=100*G[j]/float(total[j])
				percT=100*T[j]/float(total[j])
				percd=100*d[j]/float(total[j])
				perci=100*i[j]/float(total[j])
				percnc=100*nc[j]/float(total[j])
				percfs=100*fs[j]/float(total[j])
				percnon=100*non[j]/float(total[j])
				percsil=100*sil[j]/float(total[j])
				percmis=100*mis[j]/float(total[j])
				if ((j+adjustNT) % 3) == 1:
					codon_number=(j+adjustNT+2) // 3
					aa[j]=WT_AAseq[codon_number-adjustAA]+str(codon_number)
				else:
					aa[j]=""
				COUNT.write("%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n" %(j+adjustNT,WT_NTseq[j],aa[j],A[j],C[j],G[j],T[j],d[j],i[j],nc[j],total[j],percA,percC,percG,percT,percd,perci,percnc,percsil,percmis,percnon,percfs))
print("Done\n\nNucleotide mutation frequencies written to",FILENAME1,"\nNucleotide and amino acid mutations in each alignment (relative to both WT and HR GFP) written to",FILENAME2)
exit()
