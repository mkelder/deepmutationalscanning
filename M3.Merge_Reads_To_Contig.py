#!/usr/bin/python3

#quit if old python version
import sys
if sys.version_info[:2] < (3,0):
	print("This script is written for python 3.0 or above, please upgrade to a more recent version of python")
	exit()
if len(sys.argv)==2:
	INPUT = sys.argv[1]

#Set variables
ReadLen=150
AmpLen=157
#import modules
from itertools import groupby

#define functions for reverse complement and splitting numbers and text.
def complement(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)
def revcom(s):
    return complement(s[::-1])
def splitf(s):
     head = s.rstrip('0123456789')
     tail = s[len(head):]
     return head, tail

FILENAME1=INPUT+".150M.merge"
FILENAME2=INPUT+".150M.Nreads"

with open(FILENAME2,"w") as OUTFILE2: #open output FILE.
	with open(FILENAME1,"w") as OUTFILE: #open output FILE.
		with open(INPUT,"r") as SAMFILE: #open SAMFILE.
			print("Reading",INPUT)
			PREVID=""
			SEQCOUNTER=ReadCount=UnmappedCount=MappedCount=MappedPerc=UnmappedPerc=WrittenReads=WrittenN=0
			for line in SAMFILE:
					if line.startswith("@"): #skip SAM header lines	
						continue 	
					SAMfields=line.split("\t") #parse SAM fields.
					ID=(SAMfields[0]) 
					FLAG=SAMfields[1]	
					CIGAR=SAMfields[5]			
					refstart=int(SAMfields[3]) #Count where read starts wrt reference sequence.
					readNTs=list(SAMfields[9])
					readPBqual=list(SAMfields[10])
					refend=(len(readNTs)-refstart+1) #Determine length and compensate for 1-based counting
					if ("I" or "D" or "S") in CIGAR: #get rid of reads that are have indels
						continue
					#Proportions of mapped reads:
					ReadCount+=1
					if FLAG=="73" or FLAG=="137": 
						UnmappedCount+=1
					elif FLAG=="99" or FLAG=="147":
						MappedCount+=1
					if FLAG=="99": #Filter out mapped first read with mate mapped or unmapped, respectively.	
						#Scanning CIGAR string
						CIGR1=CIGAR	
						CGRS=[splitf(s) for s in CIGAR.split("S")]	
						CGRfract=[''.join(g) for _, g in groupby(CIGAR, str.isalpha)]####	
						CIGRcnt=list()
						CIGRtype=list()
						for elmt in CGRfract:
							if elmt.isdigit():
								#varno=(int(elmt)-1) 
								CIGRcnt.append(int(elmt)) #Take one to normalise to python numbering
							else:
								CIGRtype.append(elmt)		
						if CIGRcnt[0] <= 149:
							continue
						#set up empty sequence of ReadL nt, add the read nts
						NTseq=list() 
						for i in range(0,AmpLen):
							NTseq.append("-")
						rfstrt=refstart-1
						countr=RefCnt=0
						if CIGRtype[0] == 'S':
							RefCnt=CIGRcnt[0]
						CGstrt=rfstrt
						for m in CIGRtype:					
							if m == 'M':
								if CIGRcnt[countr] < 50:
									continue
								CGtemp=CGstrt	
								for n in readNTs[RefCnt:RefCnt+CIGRcnt[countr]]:
									NTseq[CGtemp]=n
									CGtemp+=1
								RefCnt=RefCnt+CIGRcnt[countr]
								CGstrt=CIGRcnt[countr]+CGstrt
							countr+=1	
						Seq1 = ''.join(str(e) for e in NTseq)
						#set up empty per base quality score of AmpLen nt
						FinalPBqual=list() 
						for i in range(0,AmpLen):
							FinalPBqual.append("!")
						rfstrt=refstart-1
						countr=RefCnt=0
						if CIGRtype[0] == 'S':
							RefCnt=CIGRcnt[0]
						CGstrt=rfstrt
						for m in CIGRtype:					
							if m == 'M': 
								CGtemp=CGstrt			
								if CIGRcnt[countr] < 30:
									continue
								for n in readPBqual[RefCnt:RefCnt+CIGRcnt[countr]]:
									FinalPBqual[CGtemp]=n
									CGtemp+=1
								RefCnt=RefCnt+CIGRcnt[countr]								
								CGstrt=CIGRcnt[countr]+CGstrt
							countr+=1		
						Qual1 = ''.join(str(e) for e in FinalPBqual)
						#Reset counters
						PREVID=ID
						SEQCOUNTER=1
					elif (ID==PREVID and FLAG=="147"): #filter out mapped second read with mate mapped or unmapped, respectively 
						if SEQCOUNTER == 2: #Check that no more than two reads are paired. 
							print("More than 2 sequences with same ID ", ID," found. Skipping read.")
							continue 	
						#Scanning CIGAR string
						CIGR2=CIGAR		
						CGRS=[splitf(s) for s in CIGAR.split("S")]			
						CGRfract=[''.join(g) for _, g in groupby(CIGAR, str.isalpha)]####	
						CIGRcnt=list()
						CIGRtype=list()
						for elmt in CGRfract:
							if elmt.isdigit():
								#varno=(int(elmt)-1) 
								CIGRcnt.append(int(elmt)) #Take one to normalise to python numbering
							else:
								CIGRtype.append(elmt)		
						if CIGRcnt[0] <= 149:
							continue
						NTseq=list() 
						for i in range(0,AmpLen):
							NTseq.append("-")
						rfstrt=refstart-1
						countr=RefCnt=0
						if CIGRtype[0] == 'S':
							RefCnt=CIGRcnt[0]
						CGstrt=rfstrt
						SeqTemp = ''.join(str(e) for e in readNTs)					
						SeqTempNTs=revcom(SeqTemp)		
						for m in CIGRtype:					
							if m == 'M': 
								CGtemp=CGstrt								
								if CIGRcnt[countr] < 50:
									continue
								for n in readNTs[RefCnt:RefCnt+CIGRcnt[countr]]:
									NTseq[CGtemp]=n
									CGtemp+=1
								RefCnt=RefCnt+CIGRcnt[countr]
								CGstrt=CIGRcnt[countr]+CGstrt
							countr+=1	
						revseq=revcom(readNTs) #not needed in dummy seqs
						Seq2 = ''.join(str(e) for e in NTseq)
						#set up empty per base quality score of AmpLen nt
						FinalPBqual=list() 
						for i in range(0,AmpLen):
							FinalPBqual.append("!")
						rfstrt=refstart-1
						countr=RefCnt=0
						if CIGRtype[0] == 'S':
							RefCnt=CIGRcnt[0]
						CGstrt=rfstrt
						for m in CIGRtype:					
							if m == 'M': 
								CGtemp=CGstrt					
								if CIGRcnt[countr] < 50:
									continue
								for n in readPBqual[RefCnt:RefCnt+CIGRcnt[countr]]:
									FinalPBqual[CGtemp]=n
									CGtemp+=1
								RefCnt=RefCnt+CIGRcnt[countr]	
								CGstrt=CIGRcnt[countr]+CGstrt
							countr+=1						
						Qual2 = ''.join(str(e) for e in FinalPBqual)
						#Compile sequences together
						if ID == 1: #In case it is the reverse, but unpaired read, it will now still output the read:	
							continue
							#OP=(ID+"    "+Seq2+"    "+Qual2+"   Unpaired"+"    "+CIGR2+"\n")
						else:
							Pbqs=list()
							SeqPB=list()
							for i in range(0,AmpLen):
								Pbqs.append("!")
								SeqPB.append("N")
							#Choose base with highest per-base quality score
							for i in range(0,AmpLen):
								if Seq1[i] == Seq2[i]:
									SeqPB[i] = Seq1[i]
									if ord(Qual1[i]) >= ord(Qual2[i]):
										Pbqs[i] = Qual1[i]
									else:
										Pbqs[i] = Qual2[i]
								else:
									Pbqs[i] = "!"
									SeqPB[i] = "N"
							for i in range(0,AmpLen):
								if Seq1[i] == Seq2[i]: #if same base is called twice, use this base with the highest quality score					
									Pbqs[i] = max(Qual1[i], Qual2[i])
									SeqPB[i] = Seq1[i]
								elif ord(Qual1[i]) >= ord(Qual2[i]) and ((ord(Qual1[i])-ord(Qual2[i])) >= 4) and ord(Qual1[i]) >= 57: #minimum difference between pbq should be 5 and higher than 57
									Pbqs[i] = Qual1[i]
									SeqPB[i] = Seq1[i].lower()
								elif ord(Qual2[i]) >= 57:
									Pbqs[i] = Qual2[i]
									SeqPB[i] = Seq2[i].lower()
								else:
									Pbqs[i] = Qual2[i]
									SeqPB[i] = "N"
							SeqFinal = ''.join(str(e) for e in SeqPB)	
							QualFinal = ''.join(str(e) for e in Pbqs)									
							OP=("@"+ID+'\n'+SeqFinal+'\n'+"+"+'\n'+QualFinal+'\n')
						if ('-') in SeqFinal:
							continue
						if ('A' or 'G' or 'C' or 'T') in SeqFinal: #check read is not empty.
							WrittenReads+=1
							if 'N' in SeqFinal[33:]:									
								OUTFILE2.write(OP)
							else:					
								WrittenN+=1
								OUTFILE.write(OP) #write to OUTFILE.
						else:
							print(SeqFinal)
						#Reset counters 
						PREVID=ID
						SEQCOUNTER=2
						ReadL1=ReadL2=0				
			#Count proportions of reads mapped
			WrittenPerc=MappedPerc=WrittenPerc=NPerc=0
			if ReadCount is not 0:
				MappedPerc=(MappedCount/ReadCount*100)
				UnmappedPerc=(UnmappedCount/ReadCount*100)
				WrittenPerc=(WrittenReads/ReadCount*100)
				NPerc=(WrittenN/ReadCount*100)
			print("Done.\n",ReadCount," reads processed.\n",MappedCount," reads mapped in pair.(",MappedPerc,"%)\n",UnmappedCount," reads mapped individually.(",UnmappedPerc,"%)\n",WrittenReads," reads merged and written.(",WrittenPerc,"%)\n",WrittenN," reads without N.(",NPerc,"%)\n")		
		SAMFILE.close()
	OUTFILE.close()
OUTFILE2.close()

