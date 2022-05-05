#!/bin/bash
### usage: ./parepare_matrix_to_velocity.sh file_path_of_expression_list process_dir script_dir
expr=$1  ##expression list from python script
process_dir=$2
script_dir=$3
expr_name=`echo $expr|awk -F "/" '{print $NF}'|awk -F "." '{print $1}'`
### step 1: filter experssion table by chosing cell number
echo "step 1: filter experssion table by chosing cell number"
  ## sort cell barcode and add the order number of cell barcode
awk '{a[$2]+=$3} END {for(i in a) print(i"\t"a[i])}' $expr | sort -k2 -n -r |awk 'BEGIN{OFS="\t"} {$3=NR;print $0}' > ${process_dir}/${expr_name}_barcode_sorted_seq.xls
Rscript ${script_dir}/cell_num_plot.R ${process_dir}/${expr_name}_barcode_sorted_seq.xls 0.05 0.05 ## Rscript to determine the number of cells.
cell_num=`cat ${process_dir}/${expr_name}_barcode_sorted_seq.xls.inflection_barcode_rank.txt` 
cut -f 1 ${process_dir}/${expr_name}_barcode_sorted_seq.xls|head -${cell_num} > ${process_dir}/cell_${cell_num}_id.txt ## cell barcode id
grep -w -Ff ${process_dir}/cell_${cell_num}_id.txt $expr > ${process_dir}/${expr_name}_${cell_num}_cell.xls ## expression list corresponding to effective cell barcode id
cut -f 1 ${process_dir}/${expr_name}_${cell_num}_cell.xls|sort|uniq > ${process_dir}/gene_${cell_num}_id.txt ## gene id
### step 2: transfer gene expression table to matrix
echo "step 2: transfer gene expression table to matrix"
awk 'ARGIND==1{a[NR]=$1} ARGIND==2{b[FNR]=$1} ARGIND==3{c[$1][$2]+=$6} END {for (i=1;i<=length(a);i++) {for (j=1;j<=length(b);j++) {printf "%d\t",c[a[i]][b[j]]};printf "\n"}}' ${process_dir}/gene_${cell_num}_id.txt ${process_dir}/cell_${cell_num}_id.txt ${process_dir}/${expr_name}_${cell_num}_cell.xls > ${process_dir}/cell_u_${cell_num}_counts.xls &  ## awk generate the matrix of unspliced counts based on cell id and gene id
awk 'ARGIND==1{a[NR]=$1} ARGIND==2{b[FNR]=$1} ARGIND==3{c[$1][$2]+=$7} END {for (i=1;i<=length(a);i++) {for (j=1;j<=length(b);j++) {printf "%d\t",c[a[i]][b[j]]};printf "\n"}}' ${process_dir}/gene_${cell_num}_id.txt ${process_dir}/cell_${cell_num}_id.txt ${process_dir}/${expr_name}_${cell_num}_cell.xls > ${process_dir}/cell_s_${cell_num}_counts.xls & ## awk generate the matrix of spliced counts based on cell id and gene id
awk 'ARGIND==1{a[NR]=$1} ARGIND==2{b[FNR]=$1} ARGIND==3{c[$1][$2]+=$4} END {for (i=1;i<=length(a);i++) {for (j=1;j<=length(b);j++) {printf "%d\t",c[a[i]][b[j]]};printf "\n"}}' ${process_dir}/gene_${cell_num}_id.txt ${process_dir}/cell_${cell_num}_id.txt ${process_dir}/${expr_name}_${cell_num}_cell.xls > ${process_dir}/cell_e_${cell_num}_counts.xls & ## awk generate the matrix of exons counts based on cell id and gene id
awk 'ARGIND==1{a[NR]=$1} ARGIND==2{b[FNR]=$1} ARGIND==3{c[$1][$2]+=$5} END {for (i=1;i<=length(a);i++) {for (j=1;j<=length(b);j++) {printf "%d\t",c[a[i]][b[j]]};printf "\n"}}' ${process_dir}/gene_${cell_num}_id.txt ${process_dir}/cell_${cell_num}_id.txt ${process_dir}/${expr_name}_${cell_num}_cell.xls > ${process_dir}/cell_i_${cell_num}_counts.xls & ## awk generate the matrix of introns counts based on cell id and gene id
awk 'ARGIND==1{a[NR]=$1} ARGIND==2{b[FNR]=$1} ARGIND==3{c[$1][$2]+=$3} END {for (i=1;i<=length(a);i++) {for (j=1;j<=length(b);j++) {printf "%d\t",c[a[i]][b[j]]};printf "\n"}}' ${process_dir}/gene_${cell_num}_id.txt ${process_dir}/cell_${cell_num}_id.txt ${process_dir}/${expr_name}_${cell_num}_cell.xls > ${process_dir}/cell_RNA_${cell_num}_counts.xls & ## awk generate the matrix of gene counts based on cell id and gene id
wait ## if the computational device did not have more than five threads, the matrix generation process should be completed one by one. 
### add col,row names to the expreesion matrix
cat ${process_dir}/cell_${cell_num}_id.txt|sed ':a;N;$!ba;s/\n/\t/g' > ${process_dir}/cell_${cell_num}_id.name
cat ${process_dir}/cell_${cell_num}_id.name ${process_dir}/cell_u_${cell_num}_counts.xls > ${process_dir}/cell_u_${cell_num}_counts.xls.col
sed "1i col/row_name" ${process_dir}/gene_${cell_num}_id.txt|paste - ${process_dir}/cell_u_${cell_num}_counts.xls.col > ${process_dir}/cell_u_${cell_num}_counts_name.xls
cat ${process_dir}/cell_${cell_num}_id.name ${process_dir}/cell_s_${cell_num}_counts.xls > ${process_dir}/cell_s_${cell_num}_counts.xls.col
sed "1i col/row_name" ${process_dir}/gene_${cell_num}_id.txt|paste - ${process_dir}/cell_s_${cell_num}_counts.xls.col > ${process_dir}/cell_s_${cell_num}_counts_name.xls
cat ${process_dir}/cell_${cell_num}_id.name ${process_dir}/cell_e_${cell_num}_counts.xls > ${process_dir}/cell_e_${cell_num}_counts.xls.col
sed "1i col/row_name" ${process_dir}/gene_${cell_num}_id.txt|paste - ${process_dir}/cell_e_${cell_num}_counts.xls.col > ${process_dir}/cell_e_${cell_num}_counts_name.xls
cat ${process_dir}/cell_${cell_num}_id.name ${process_dir}/cell_i_${cell_num}_counts.xls > ${process_dir}/cell_i_${cell_num}_counts.xls.col
sed "1i col/row_name" ${process_dir}/gene_${cell_num}_id.txt|paste - ${process_dir}/cell_i_${cell_num}_counts.xls.col > ${process_dir}/cell_i_${cell_num}_counts_name.xls
cat ${process_dir}/cell_${cell_num}_id.name ${process_dir}/cell_RNA_${cell_num}_counts.xls > ${process_dir}/cell_RNA_${cell_num}_counts.xls.col
sed "1i col/row_name" ${process_dir}/gene_${cell_num}_id.txt|paste - ${process_dir}/cell_RNA_${cell_num}_counts.xls.col > ${process_dir}/cell_RNA_${cell_num}_counts_name.xls
rm ${process_dir}/*.col ${process_dir}/cell_${cell_num}_id.name ${process_dir}/cell_*_${cell_num}_counts.xls
