name=b48h1-mht

run_path=/project/$name

db="/database/Macaca_fascicularis/NCBI/GCF_000364345.1/all.fa"

/zfswh1/BC_RD_P0/heyingdong/project/IsoSeq_SingelCell/ProjectMouseSperm/software/minimap2/minimap2 -ax splice -t 16 -uf --secondary=no -C5 $db  ${run_path}/FL/fl_seq.fq > ${run_path}/FL/fl_seq.sam
