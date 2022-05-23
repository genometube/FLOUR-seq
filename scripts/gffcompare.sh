#git clone https://github.com/gpertea/gscripts.git
name=b48h1-mht
run_path=/project/$name/process/

perl /software/gscripts-master/sam2gff.pl ${run_path}/FL/fl_seq.sam |grep "\.m1"|sed 's#\.m[0-9]##g' > ${run_path}/FL/gene/gffcompare/all.mapping.gff
/software/gffcompare-0.11.6/gffcompare -r /database/GTF/mRNA.gtf  ${run_path}/FL/gene/gffcompare/all.mapping.gff  -o ${run_path}/FL/gene/gffcompare/gffcompare


