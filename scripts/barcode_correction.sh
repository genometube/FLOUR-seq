name=b48h1-mht
run_path=/project/$name
mkdir -p /project/$name/bc_correct
cd /project/$name/bc_correct
## barcode correction
sed 'N;N;N;s#\n#\t#g'  ../FL/isoseq_flnc.BarcodeUMI.fastq|awk '{bc=substr($2,1,16);umi=substr($2,17,28);bc_qua=substr($4,1,16);umi_qua=substr($2,17,28);print $1"\n"bc"\n+\n"bc_qua}' > /project/$name/FL/Barcode.fastq &

sed 'N;N;N;s#\n#\t#g' ../FL/isoseq_flnc.BarcodeUMI.fastq|awk '{bc=substr($2,1,16);umi=substr($2,17,28);bc_qua=substr($4,1,16);umi_qua=substr($2,17,28);print $1"\n"umi"\n+\n"umi_qua}' > /project/$name/FL/umi.fa


/software/cellranger-3.1.0/miniconda-cr-cs/4.3.21-miniconda-cr-cs-c10/bin/python  /script/bc_corrector.py  -w  /barcode/3M-february-2018.txt -t 0.95 -i ${run_path}/FL/Barcode.fastq -o  /project/$name/FL/