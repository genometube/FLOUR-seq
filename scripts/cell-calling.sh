script_path="/script"
software_path="/software"
name=b48h1-mht
run_path=/project/$name

process_dir="${run_path}/process/"
shell_dir="${run_path}/shell"

mkdir -p ${process_dir}/FL/gene/gffcompare
# mkdir -p ${process_dir}/FL/gene/gffcompare

cd  ${shell_dir}/FL/process

#gene_bc.mtx format
#ENSMUSG00000003429	AAAAAAAAAAAAAAAA	1

#--genome mm10 / hg19 / ...

awk 'NR>1' ${process_dir}/FL/gene/gffcompare/gene_expr.xls|cut -f1-3 > ${process_dir}/FL/gene/gene_bc.mtx
${software_path}/Rscript ${script_path}/callcell_cr.R -i ${process_dir}/FL/gene/gene_bc.mtx -o ${process_dir}/FL/gene/ -e 10000 --genome monkey

