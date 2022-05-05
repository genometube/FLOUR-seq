#!/bin/bash
#dependencies: fastp, seqtk, bowtie2, samtools, python3
#Example: for test.fastq.gz, run "sh debarcode_HISOFA.sh test"
data="$1.fastq.gz"
mkdir -p $1_out/tmp
cd $1_out 

##s0_cutP5P7
fastp -i ${data} -o $1_1.fastq.gz --adapter_sequence ATGTACTCTGCGTTG --disable_quality_filtering &&\
fastp -i $1_1.fastq.gz -o $1_2.fastq.gz --adapter_sequence AGGGAAAGAGTGTAG --disable_quality_filtering &&\
seqtk seq -r $1_2.fastq.gz |gzip > $1_3.fastq.gz &&\
rm $1_1.fastq.gz $1_2.fastq.gz &&\
fastp -i $1_3.fastq.gz -o $1_4.fastq.gz --adapter_sequence ATGTACTCTGCGTTG --disable_quality_filtering &&\
fastp -i $1_4.fastq.gz -o $1_5.fastq.gz --adapter_sequence AGGGAAAGAGTGTAG --disable_quality_filtering &&\
seqtk seq -r $1_5.fastq.gz |gzip > $1_6.fastq.gz &&\
rm $1_3.fastq.gz $1_4.fastq.gz $1_5.fastq.gz fastp.html fastp.json &&\

##s0_SE300_reads
zcat $1_6.fastq.gz |perl -nlE 'if((length $_)>300){say substr($_,0,150).substr($_,-150,150)} else{say $_}' |gzip > "SE300_$1.fastq.gz" &&\

##s0_SE300_sam
bowtie2 -x cell_label.fa -U "SE300_$1.fastq.gz" -S $1.sam --local --np 0 --mp 1,0 --rdg 0,1 --rfg 0,1 --score-min L,152,0 -L 23 -D 200 -R 30 -i S,1,0.2 2>bowtie2.log &&\
rm "SE300_$1.fastq.gz" &&\

##s1_SE300_forward & s1_SE300_reverse
samtools view -F4 -f16 $1.sam |cut -f1,3 > "reverse_$1.csv" &&\
samtools view -F4 -F16 $1.sam |cut -f1,3 > "forward_$1.csv" &&\
rm $1.sam &&\

##s2_SE300_for_fq & s2_SE300_rev_fq
cut -f1 "forward_$1.csv" > "tmp/fwd_ID_$1.csv" &&\
seqtk subseq "$1_6.fastq.gz" "tmp/fwd_ID_$1.csv" > "fwd_fq_$1.fq" &&\
cut -f1 "reverse_$1.csv" > "tmp/rev_ID_$1.csv" &&\
seqtk subseq "$1_6.fastq.gz" "tmp/rev_ID_$1.csv" > "rev_fq_$1.fq" &&\
rm "tmp/fwd_ID_$1.csv" "tmp/rev_ID_$1.csv" "$1_6.fastq.gz" &&\

##s3_SE300_for_pT
python Find_polyA.py --fastq "fwd_fq_$1.fq" --polyA no --polyT yes > "tmp/for_pT_s1_$1.csv" &&\
awk '{print $1,$(NF)}' "tmp/for_pT_s1_$1.csv" |awk -F[\s,] '{print $1,$2}' |awk '{print $1,$3}' > "for_pT_$1.csv" &&\
##s3_SE300_rev_pA
python Find_polyA.py --fastq "rev_fq_$1.fq" --polyA yes --polyT no > "tmp/rev_pT_s1_$1.csv" &&\
awk '{print $1"\t"$(NF)}' "tmp/rev_pT_s1_$1.csv" |sed -r 's/;/\t/g' |awk '{print $1,$NF}' |awk -F, '{print $1}' > "rev_pA_$1.csv" &&\
rm "tmp/for_pT_s1_$1.csv" "tmp/rev_pT_s1_$1.csv" &&\

##s4_SE300_cut_pA & s4_SE300_cut_pT
python cut_pT_pA.py "rev_pA_$1.csv" "rev_fq_$1.fq" "cutpA_$1.fq" polyA &&\
python cut_pT_pA.py "for_pT_$1.csv" "fwd_fq_$1.fq" "cutpT_$1.fq" polyT &&\
gzip "cutpA_$1.fq" &&\
gzip "cutpT_$1.fq" &&\
rm "for_pT_$1.csv" "rev_pA_$1.csv" "fwd_fq_$1.fq" "rev_fq_$1.fq" "noA.fq" "noT.fq"&&\
echo ========== Lib "$1" is done ==========
