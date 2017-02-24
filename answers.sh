#! /usr/bin/env bash


#1: Use BEDtools intersect to identify the size of the largest overlap 
# between CTCF and H3K4me3 locations

CTCF="$HOME/Genome_Class/rc-files/data-sets/bed/encode.tfbs.chr22.bed.gz"

gzcat $CTCF | awk '($4=="CTCF") {print $0}' > CTCF.bed

H3K="$HOME/Genome_Class/rc-files/data-sets/bed/encode.h3k4me3.hela.chr22.bed.gz"

answer_1=$(bedtools intersect -a CTCF.bed -b $H3K \
    | awk 'BEGIN{OFS="\t"} {print $0, $3-$2}' \
    | sort -k5nr \
    | cut -f5 \
    | head -n1)

echo "answer-1: $answer_1"



#2: Use BEDtools to calculate the GC content of nucleotides 19,000,000 to
# 19,000,500 on chr22 of hg19 genome build. Report the GC content as a fraction.

fastseq="$HOME/Genome_Class/rc-files/data-sets/fasta/hg19.chr22.fa"

answer_2=$(echo -e "chr22\t19000000\t19000500" > interval.bed \
 | bedtools nuc -fi $fastseq -bed interval.bed \
 | awk 'BEGIN {OFS="\t"} {print $6, $7, $8, $9}' \
 | tail -n 1 \
 | awk '{print ($2+$3)/($1+$2+$3+$4)}')

echo "answer-2: $answer_2"




#3 Use BEDtools to identify the length of the CTCF ChIP-seq peak (i.e., interval) 
# that has the largest mean signal in ctcf.hela.chr22.bg.gz.

Bedgraph="$HOME/Genome_Class/rc-files/data-sets/bedtools/ctcf.hela.chr22.bg.gz"

answer_3=$(bedtools map -c 4 -o mean -a CTCF.bed -b $Bedgraph \
    | sort -k5nr \
    | awk 'BEGIN {OFS="\t"} {print $0, $3-$2}' \
    | head -n 1 \
    | awk '{print $6}')

echo "answer-3: $answer_3"




#4 Use BEDtools to identify the gene promoter (defined as 1000 bp upstream
# of a TSS) with the highest median signal in ctcf.hela.chr22.bg.gz. Report the gene name 

TSS="$HOME/Genome_Class/rc-files/data-sets/bed/tss.hg19.chr22.bed"

Genome="$HOME/Genome_Class/rc-files/data-sets/genome/hg19.genome"

CTCF_4="$HOME/Genome_Class/rc-files/data-sets/bedtools/ctcf.hela.chr22.bg.gz"

answer_4=$(bedtools flank -i $TSS -g $Genome -l 1000 -r 0 > flanked.bed 
    bedtools sort -i flanked.bed > sorted.bed 
    bedtools map -c 4 -o median -a sorted.bed -b $CTCF_4 \
    | sort -k7nr \
    | cut -f 4 \
    | head -n 1)

echo "answer-4: $answer_4"    



#5 Use BEDtools to identify the longest interval on chr22 that is not
# covered by genes.hg19.bed.gz. Report the interval like chr1:100-500

hg19genes="$HOME/Genome_Class/rc-files/data-sets/bed/genes.hg19.bed"

answer_5=$(bedtools complement -i $hg19genes -g $Genome \
    | awk '($1="chr22") {print $0, $3-$2}' \
    | sort -k4nr \
    | head -n 1 \
    | awk '{print $1 ":" $2 "-" $3}')

echo "answer-5: $answer_5"



#6 EC. Use bedtools genomcov to see the coverage of h3k4me
# sites on chr22 throughout hg19.genome

# H3="$HOME/Genome_Class/rc-files/data-sets/bed/encode.h3k4me3.hela.chr22.bed.gz"

# answer_6=$(bedtools sort -i $H3 > sortedH3.bed
 #   bedtools genomecov -i sortedH3.bed -g $hg19genes -bg \
  #  | sort -k4nr \
   # | head -n 20 )

# echo "answer-6: $answer_6"















