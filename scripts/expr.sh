script_path="/script"
software_path="/software"
name=b48h1-mht
run_path=/project/$name

barcode=/project/$name/FL/Barcode.fastq
umi=/project/$name/FL/umi.fa

process_dir="${run_path}/process"
shell_dir="${run_path}/shell"

mkdir -p ${shell_dir}/FL/gene/

mkdir -p ${process_dir}/FL/gene/gffcompare
# mkdir -p ${process_dir}/FL/gene/gffcompare

cd  ${shell_dir}/FL/gene

##= c k m n j e o I y
awk 'NR==FNR{a[$1]=$2"\t"$3;next}{if($4 in a)print $0"\t"a[$4]}'  ${process_dir}/exon_tron/gene/exon_intron.stat_latest.xls ${process_dir}/FL/gene/gffcompare/gffcompare.all.mapping.gff.tmap| awk '$3=="=" || $3=="u" || $3=="k"|| $3=="m"|| $3=="n"|| $3=="j"|| $3=="e"|| $3=="o"|| $3=="i"|| $3=="y" ' > ${process_dir}/FL/gene/gffcompare/gffcompare.all.mapping.gff.tmap.2


/software/cellranger-3.1.0/miniconda-cr-cs/4.3.21-miniconda-cr-cs-c10/bin/python /script/gene_BC_mtx2_exon_intron.py --input ${process_dir}/FL/gene/gffcompare/gffcompare.all.mapping.gff.tmap.2 --barcode $barcode --umi $umi > ${process_dir}/FL/gene/gffcompare/gene_expr.xls 

