software="/software/"
script_path="/script"
name=b48h1-mht
##################

## mm10


db="/database/Macaca_fascicularis/NCBI/GCF_000364345.1/"
#GTF/mRNA.gtf


run_path="/zfswh1/BC_RD_P0/P19Z12200N0089/2021/BD_ONT/sc_isoseq/run/run_2022_03/b48h1-mht/"
process_dir="${run_path}/process/$name"
shell_dir="${run_path}/shell/$name"

##extract exon _intron
process_dir="${run_path}/process/exon_tron/gene"
mkdir -p ${process_dir}

# awk '$3=="exon"' $db/genome/refdata-cellranger-mm10-3.0.0/genes/genes.gtf  > $db/refdata-cellranger-mm10-3.0.0/genes/exon_genes.gtf

awk '$3=="exon"' $db/GTF/mRNA.gtf |sort -k1,1 -k4,4n    > ${process_dir}/exon_genes.gtf

python ${script_path}/exon_gtf.py ${process_dir}/exon_genes.gtf |sort -k1,1 -k4,4n|awk '{print $1"\t"$4"\t"$5"\t"$2}'  > ${process_dir}/exon.bed

#merge exon
#samtools faidx chrALL.fa > chrALL.fa.fai
cat $db/GenomeGatkIndex/chrALL.fa.fai|sort -k1,1 -k2,2n > ${process_dir}/chrALL.fa.fai

awk '$4=="exon_DIY"' ${process_dir}/exon.bed|$software/bedtools merge -i - |awk '{print $1"\t"$2"\t"$3"\texon_DIY"}' |sort -k1,1 -k2,2n > ${process_dir}/exon_merge.bed
# awk '$4=="exon_DIY"' ${process_dir}/exon.bed|$software/bedtools merge -i - |awk '{print $1"\t"$2"\t"$3"\texon_DIY"}' |sort -k1,1 -k2,2n|grep "chr"> ${process_dir}/exon_merge.bed

$software/bedtools complement -i ${process_dir}/exon_merge.bed  -g  ${process_dir}/chrALL.fa.fai |awk '{print $0"\tintron_DIY"}'>  ${process_dir}/intron.bed

cat ${process_dir}/exon.bed ${process_dir}/intron.bed|sort -k1,1 -k2,2n |uniq > ${process_dir}/merge_exon_intron.uniq.bed

awk '$3=="exon"' ${process_dir}/../../FL/gene/gffcompare/all.mapping.gff > ${process_dir}/../../FL/gene/all.collapsed.gff2

awk '{print $1"\t"$4"\t"$5"\t"$9}' ${process_dir}/../../FL/gene/all.collapsed.gff2 |sed 's#Parent=##g;s#\.m[0-9]##g' |sort|uniq |sort -k1,1 -k3,3n > ${process_dir}/fl_seq.sorted.bed


$software/bedtools intersect -a ${process_dir}/merge_exon_intron.uniq.bed  -b ${process_dir}/fl_seq.sorted.bed  -wo  > ${process_dir}/intersect_bed

#get longest uniq intersect 
sort -k5,5 -k8,8 -k6,6n -k9,9nr ${process_dir}/intersect_bed |awk '!($4$6$7$8 in a){a[$4$6$7$8];print $0}' >  ${process_dir}/intersect_bed.uniq 

echo -e "gene\texon\tintron" > ${process_dir}/exon_intron.stat_latest.xls;

python ${script_path}/extract_intron_exon.py  ${process_dir}/intersect_bed.uniq  |sed 's#"##g;s#;##g'>>  ${process_dir}/exon_intron.stat_latest.xls