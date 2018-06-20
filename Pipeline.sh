#!/bin/bash -f
#Grid Engine options: 4gb in 8 cores, expected run time less than 10 hours
#$ -l h_vmem=4g
#$ -l h_rt=10:00:00
#$ -pe sharedmem 8
#check amount of memory seen by jobqstat
ulimit -v
#Configure modules
. /etc/profile.d/modules.sh

#Prior to the analyses, several reference files need to be in place. 
#Place the FASTQ files compressed in the GNU zip (indicated by the .gz file extension) in analysis directory.
#Place the reference FASTA file, containing the sequence name in the first line in format “>SEQUENCE NAME” and the reference sequence in the second line, in the analysis directory.
#Create a filenames_generic.txt file containing the sample names in the analysis directory.
###Set parameters below prior to running file
DIR="/mnt/nfsstaging/igmm/datastore/awood-lab/Martijn/New/"
FASTA="GFP_WT.fasta"
HRCORE="SILENT" #SILENT or STOP

###Turn off/on modules 6-12 by toggling 0/1. Modules 1-5 are essential and on by default.
MOD6=1
MOD7=1
MOD8=1
MOD9=1
MOD10=1
MOD11=1
MOD12=1

#Check parameters are valid and navigate to directory.
	if [ -z "$DIR" ]; then 
		echo "No analysis directory submitted. exiting script." 
		exit; else 
		cd "$DIR"
	fi
	if [ -f $FASTA ]; then
		echo "FASTA file found, continuing. "
	else
		echo "Target directory must contain the indicated FASTA file,  but file was not found. Exiting script."
	fi

#Load modules
module add igmm/apps/bowtie/2.3.1
module add igmm/apps/TrimGalore/0.4.1
module add igmm/apps/python/3.5.2
module add igmm/apps/R/3.3.3 
module add igmm/apps/weblogo/3.5.0

#Unzippinf all gzips
for a in *.gz
do
  gunzip $a
done

#Store consensus sequence in variable
CONSENSUS=$(cat $FASTA | head -2 | tail -1)

####Module 1: Trimming reads (Essential)
#This module will run Trim Galore! (a wrapper for cutadapt and fastqc) on the list of files in filenames.txt in the target directory.
#filenames.txt will automatically be generated based on the .fastq files in the directory.
#In paired end mode the first two files in filenames.txt are paired, then the next two etc.
	#get filenames and run trimgalore
	echo "####Module 1: Trimming Reads (Essential)"
	ls -1 *.fastq > filenames.txt
	filecount=($(wc filenames.txt)[0])
	i=1
	j=$[$filecount/2]

	if [ -f filenames.txt ]; then
			{ IFS=$'\n'; FILENAMES=($(cat filenames.txt)); }
	else
		echo "Target directory must contain the file filenames.txt containing a list of files to process. Exiting script."
		exit
	fi

	while [ "$i" -le "$j" ]; do
		FILEINDEX=$[2*$i]
		FILEINDEX=$[$FILEINDEX-1]
		FILE2=${FILENAMES["$FILEINDEX"]}
		FILEINDEX=$[$FILEINDEX-1]
		FILE1=${FILENAMES["$FILEINDEX"]}
					trim_galore --stringency 4 --paired --nextera "$FILE1" "$FILE2"
		FILEROOT1=${FILE1%.*}
		if [ -f "$FILEROOT1".sanfastq_val_1.fq   ]; then mv "$FILEROOT1".sanfastq_val_1.fq cut."$FILE1"; fi
		FILEROOT2=${FILE2%.*}
		if [ -f "$FILEROOT2".sanfastq_val_2.fq ]; then mv "$FILEROOT2".sanfastq_val_2.fq cut."$FILE2"; fi
		i=$[$i+1]
	done
	#Grep mapping efficiencies, saves statistics to respective txt files.
	mkdir analysis
	grep 'Total written' *.txt > analysis/Totalwritten.txt
	grep 'Total reads processed' *.txt > analysis/Totalprocessed.txt
	grep 'Total written (filtered)' *.txt > analysis/Totalfiltered.txt
	echo "Finished with module 1."
	echo "-------------------------------------------------------------------------"

####Module 2: Aligning Reads to Reference 	(Essential)
#Module will run BowTie2 on the list of files in filenames_trimmed.txt in the target directory using the fasta file set at the start of the script.
#filenames_trimmed.txt will automatically be generated based on the .fq files in the directory.
#Mismatch penalty (mp) can be increased to make substitutions more stringent
	echo "####Module 2: Aligning Reads to Reference 	(Essential)"
	#Build Bowtie2 index
	echo "Building index"
	bowtie2-build "scripts/$FASTA" "$FASTA"

	#Get filenames and run Bowtie2
	ls -1 *.fq > filenames_trimmed.txt
	i=1

	if [ -f filenames_trimmed.txt ]; then
		{ IFS=$'\n'; FILENAMES=($(cat filenames_trimmed.txt)); }
	else
		echo "Target directory must contain the file filenames_trimmed.txt containing a list of files to process. Exiting script."
		exit
	fi

	while [ "$i" -le "$j" ]; do
		FILEINDEX=$[2*$i]
		FILEINDEX=$[$FILEINDEX-1]
		FILE2=${FILENAMES["$FILEINDEX"]?_val_2.fq}
		FILEINDEX=$[$FILEINDEX-1]
		FILE1=${FILENAMES["$FILEINDEX"]?_val_1.fq}
					echo "Mate 1: $FILE1"
					echo "Mate 2: $FILE2"
					bowtie2 --threads 8 --phred33  --no-mixed --mp 2 -i G,1,0.25 --rdg 9,3 --rfg 9,3 -N 1 --rg-id 3 -x GFP_WT.fasta -1 $FILE1 -2 $FILE2 -S $FILE1.sam
		FILEROOT1=${FILE1%.*}
		if [ -f "$FILEROOT1".sanfastq_val_1.fq   ]; then mv "$FILEROOT1".sanfastq_val_1.fq cut."$FILE1"; fi
		FILEROOT2=${FILE2%.*}
		if [ -f "$FILEROOT2".sanfastq_val_2.fq ]; then mv "$FILEROOT2".sanfastq_val_2.fq cut."$FILE2"; fi
		i=$[$i+1]
	done
	echo "Finished with module 2."
	echo "-------------------------------------------------------------------------"

####Module 3: Merging reads into single contig 	(Essential)
#This module merges forward and reverse reads into a single contig of predefined amplicon length.
#Adjust parameters in 'Merge_Reads_To_Contig.py' to adjust read length, amplicon length, stringencies and the toggle the inclusion of indels and ambiguous bases
	echo "####Module 3: Merging reads into single contig"
	for file in *.sam; do
		python3 scripts/M3.Merge_Reads_To_Contig.py $file
	done
	echo "Finished with module 3."
	echo "-------------------------------------------------------------------------"

####Module 4: Remapping reads into SAM files 	(Essential)
#This module remaps the merged contigs to the reference sequence so that the files are in the SAM format again. 
	echo "####Module 4: Remapping reads into SAM files"
	for file in *.merge; do
		FILEROOT=${file%.*}
		bowtie2 -x GFP_WT.fasta -U $file -S $FILEROOT.remapped.sam
		echo "saved $FILEROOT"
	done
	echo "Finished with module 4."
	echo "-------------------------------------------------------------------------"

####Module 5: Sorting sequences based on core sequence 	(Essential)	
#This module sorts sequences based on their core sequence into three classes: NOHR (i.e. wildtype), HR and AMB (ambiguous)
#Core frequencies will be displayed on the command line
#Reference sequence and core positions can be set at the start of the script
	echo "####Module 5: Sorting sequences based on core sequence"
	for file in *remapped.sam; do
		python3 scripts/M5.Sort_SAM.py $file $HRCORE $CONSENSUS
	done
	echo "Finished with module 5."
	echo "-------------------------------------------------------------------------"

#Move unused files into respective directories and only retain reads with HR cores
	echo "Sorting unused files into respective folders"
	mkdir fastq trimgalore_reports trimmed samfiles1 samfiles_remapped merge_Nreads merge_reads NOHR AMB
	mv *.fastq fastq
	mv *_trimming_report.txt trimgalore_reports
	mv *.fq trimmed
	mv *remapped.sam samfiles_remapped
	mv *.sam samfiles_remapped
	mv *.Nreads merge_Nreads
	mv *.merge merge_reads
	mv *.NOHR NOHR
	mv *.AMB AMB
	echo "Finished file sorting."
	echo "-------------------------------------------------------------------------"

####Module 6: Check mutation frequencies in WT reads (Optional)	
#This script scans wildtype core and assesses the proportion of wildtype reads with full consensus. 
if [ $MOD6 == 1 ];then
	echo "####Module 6: Check mutation frequencies in WT reads"
	cd NOHR
	for file in *.NOHR
	do
	python3 ../scripts/M6.Filter_WT_mutations.py $file $CONSENSUS
	done
	cd ..
	echo "Finished with module 6."
	echo "-------------------------------------------------------------------------"
fi
	
####Module 7: Filter HR reads based on number of nucleotide changes (Optional)
#Filter HR reads based on the number of non-consensus nucleotides, discard all other HR reads
if [ $MOD7 == 1 ];then
	echo "####Module 7: Filter HR reads based on number of nucleotide changes"
	for file in *$core.HR; do
		python3 scripts/M7.Filter_HR_uniqueNT.py $file $HRCORE $CONSENSUS
	done
	echo "Finished with module 7."
	echo "-------------------------------------------------------------------------"
fi

####Module 8: Counting nucleotide frequencies and mutations 	(Optional)
#This module scans through the (sorted) SAM file and counts nucleotide frequencies for each base position across all the reads and outputs these into *.NTcount
#It furthermore converts nucleotide variants into amino acid substitutions and outputs all reads with their amino acid substitutions into *.mutations
#Reference sequence and core positions can be set at the start of the script
if [ $MOD8 == 1 ];then
	echo "####Module 8: Counting nucleotide frequencies and mutations"
	for file in *$core.HR; do #Can alternatively be adjusted to HR/NOHR/AMB
		python3 scripts/M8.Mutation_SAM.py $file $HRCORE $CONSENSUS
	done
	echo "Finished with module 8."
	echo "-------------------------------------------------------------------------"
fi

####Module 9: Codon count  	(Optional)	
#This module counts each of the codons in the predefined reading frame
#Output is a table of frequency of each possible codon at each position
#Reading frame and start amino acid need to be preset at the beginning of the script
if [ $MOD9 == 1 ];then
	echo "####Module 9: Plotting sequence logos"
	for file in *.HR; do
		python3 scripts/M9.Codoncount.py $file $HRCORE $CONSENSUS
	done
	echo "Finished with module 9."
	echo "-------------------------------------------------------------------------"

fi

####Module 10: Plotting non-consensus nucleotide frequencies  	(Optional)	
#The module creates non-consensus nucleotide scatter plots based on the nucleotide frequencies in *.NTcount
#The additional script plots stacked barplot of nucleotide substituion frequencies (specific for GFP project)
#Optimised for GFP, positional parameters can be adjusted to use it for other positions. 
if [ $MOD10 == 1 ];then
	echo "####Module 10: Plotting non-consensus nucleotide frequencies"
	ls -1 *.fastq > filenames_generic.txt
	FILENAMES=($(cat filenames_generic.txt))
	for class in HR; do #Can alternatively be adjusted to exclude any of HR/NOHR/AMB
	i=0
	for file in *$HRCORE.$class.NTcount; do
	FILEROOT=${file%.*}
	FNAME=${FILENAMES["$i"]}
	echo $FILEROOT
	echo $FNAME
	Rscript scripts/M10.Plot_nucleotide_scatter.R $FNAME $file
	Rscript scripts/M10b.GFP_stacked_barplot.R $FNAME $file
	echo " "
	i=$[$i+1]
	done
	done
	#fix cryptic file names
	for file in *; do mv "$file" "$(echo $file | sed s'/\r//g')"; done
	echo "Finished with module 10."
	echo "-------------------------------------------------------------------------"

fi

####Module 11: Count frequency of single amino acid changes (Optional)
#This module requires 'GFP_ddG_values.csv' to be present in the directory
#This module takes HR reads with only a single amino acid change and outputs a table of amino acid changes per sample
#It also outputs a table with counts for amino acid changes in the variable positions only (".var")
#Looks up ddG values calculated for GFP in CSV table (specific for GFP project)
#Gives an error in case the *.mutations file is empty
if [ $MOD11 == 1 ];then
	echo "####Module 11: Count frequency of single amino acid changes"
	for file in *.mutations; do
	echo $file
	Rscript scripts/M11.Sort_uniqueAAmuts.R $file
	echo " "
	done
	mv *.CSV analysis
	echo "Finished with module 11."
	echo "-------------------------------------------------------------------------"

fi

####Module 12: Plotting sequence logos (Optional)
#This module converts HR contigs to amino acid sequences and plotting sequence logos
if [ $MOD12 == 1 ];then
	echo "####Module 12: Plotting sequence logos"
	for file in *.HR; do
	python3 scripts/M12.SAMtoAA.py $file $HRCORE $file $file $file
	done
	module add igmm/apps/weblogo/3.5.0
	for file in *.AAs; do
	echo $file
	weblogo -F pdf -n  80 -i 83 --errorbars NO -l 97 -u 126 -s large -c chemistry < $file > $file.pdf
	done 
	mkdir AA
	mv *.AAs AA
	mv *.pdf analysis
	echo "Finished with module 12."
fi
echo ""
echo "-------------------------------------------------------------------------"
echo "Finished with analyses. Output files can be found in /analysis directory."
echo "-------------------------------------------------------------------------"
echo ""
