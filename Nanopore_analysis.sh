name="sperm"
software="/software"
run_path=/project/$name ##The path the result from debarcode_HISOFA.sh
db="/database/Macaca_fascicularis/NCBI/GCF_000364345.1/"

### Mapping
db="/database/Macaca_fascicularis/NCBI/GCF_000364345.1/all.fa"
zcat ${run_path}/*fq.gz > ${run_path}/fl_seq.fq
zcat ${run_path}/*.csv.gz > ${run_path}/barcode.xls
${software}/minimap2 -ax splice -t 16 -uf --secondary=no -C5 $db  ${run_path}/fl_seq.fq > ${run_path}/fl_seq.sam

### Mapping reads to reference annotation
#git clone https://github.com/gpertea/gscripts.git to find required script
mkdir -p ${run_path}/process

perl /software/gscripts-master/sam2gff.pl ${run_path}/fl_seq.sam |grep "\.m1"|sed 's#\.m[0-9]##g' > ${run_path}/process/gffcompare/all.mapping.gff
/software/gffcompare-0.11.6/gffcompare -r $db/GTF/mRNA.gtf  ${run_path}/process/gffcompare/all.mapping.gff  -o ${run_path}/process/gffcompare/gffcompare

### Extract exon intron
process_dir="${run_path}/process"
mkdir -p $process_dir
awk '$3=="exon"' $db/GTF/mRNA.gtf |sort -k1,1 -k4,4n    > ${process_dir}/exon_genes.gtf
python ${script_path}/exon_gtf.py ${process_dir}/exon_genes.gtf |sort -k1,1 -k4,4n|awk '{print $1"\t"$4"\t"$5"\t"$2}'  > ${process_dir}/exon.bed
## Merge exon
#samtools faidx chrALL.fa > chrALL.fa.fai
cat $db/GenomeGatkIndex/chrALL.fa.fai|sort -k1,1 -k2,2n > ${process_dir}/chrALL.fa.fai
awk '$4=="exon_DIY"' ${process_dir}/exon.bed|$software/bedtools merge -i - |awk '{print $1"\t"$2"\t"$3"\texon_DIY"}' |sort -k1,1 -k2,2n > ${process_dir}/exon_merge.bed
$software/bedtools complement -i ${process_dir}/exon_merge.bed  -g  ${process_dir}/chrALL.fa.fai |awk '{print $0"\tintron_DIY"}'>  ${process_dir}/intron.bed
cat ${process_dir}/exon.bed ${process_dir}/intron.bed|sort -k1,1 -k2,2n |uniq > ${process_dir}/merge_exon_intron.uniq.bed
awk '$3=="exon"' ${process_dir}/gffcompare/all.mapping.gff > ${process_dir}/all.collapsed.gff2
awk '{print $1"\t"$4"\t"$5"\t"$9}' ${process_dir}/all.collapsed.gff2 |sed 's#Parent=##g;s#\.m[0-9]##g' |sort|uniq |sort -k1,1 -k3,3n > ${process_dir}/fl_seq.sorted.bed
$software/bedtools intersect -a ${process_dir}/merge_exon_intron.uniq.bed  -b ${process_dir}/fl_seq.sorted.bed  -wo  > ${process_dir}/intersect_bed
## Get longest uniq intersect
sort -k5,5 -k8,8 -k6,6n -k9,9nr ${process_dir}/intersect_bed |awk '!($4$6$7$8 in a){a[$4$6$7$8];print $0}' >  ${process_dir}/intersect_bed.uniq
echo -e "gene\texon\tintron" > ${process_dir}/exon_intron.stat_latest.xls
python ${script_path}/extract_intron_exon.py  ${process_dir}/intersect_bed.uniq  |sed 's#"##g;s#;##g'>>  ${process_dir}/exon_intron.stat_latest.xls
### Generate gene-barcode exon intron-cell matrix
##= c k m n j e o I y
awk 'NR==FNR{a[$1]=$2"\t"$3;next}{if($4 in a)print $0"\t"a[$4]}'  ${process_dir}/exon_intron.stat_latest.xls ${process_dir}/gffcompare/gffcompare.all.mapping.gff.tmap| awk '$3=="=" || $3=="u" || $3=="k"|| $3=="m"|| $3=="n"|| $3=="j"|| $3=="e"|| $3=="o"|| $3=="i"|| $3=="y" ' > ${process_dir}/gffcompare/gffcompare.all.mapping.gff.tmap.2
awk '{print ">"$1"\n"$2"\n""+""\n""+"}' ${run_path}/barcode.xls > ${run_path}/barcode.fq
awk '{print ">"$1"\n"NR}' ${run_path}/barcode.xls > ${run_path}/umi.fa
python ${script_path}/gene_BC_mtx2_exon_intron.py --input ${process_dir}/gffcompare/gffcompare.all.mapping.gff.tmap.2 --barcode ${run_path}/barcode.fq --umi ${run_path}/umi.fa > ${process_dir}/gffcompare/gene_expr.xls 
### Cell-calling using gene-barcode matrix
# gene_bc.mtx format
# ENSMUSG00000003429	AAAAAAAAAAAAAAAA	1
awk 'NR>1' ${process_dir}/gffcompare/gene_expr.xls|cut -f1-3 > ${process_dir}/gene_bc.mtx
${software}/Rscript ${script_path}/callcell_cr.R -i ${process_dir}/gene_bc.mtx -o ${process_dir}/ -e 10000 --genome mouse
