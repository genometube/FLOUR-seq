name=b48h1-mht
bam=/project/RNA_velocity/HIT_data/data/${name}-ccs.bam

primer_db=/project/IsoSeq_SingelCell/ProjectMouseSperm/db/Primer.fa
run_path=/project/$name

db="/database/Macaca_fascicularis/NCBI/GCF_000364345.1/all.fa"

export LD_LIBRARY_PATH=/share/app/glibc-2.14/lib:$LD_LIBRARY_PATH
source /ifswh1/BC_PUB/biosoft/BC_NQ/01.Soft/environment.sh

## download scISA-Tools
##git clone https://github.com/shizhuoxing/scISA-Tools.git
script_path=/software/scISA-Tools/bin/



mkdir -p ${run_path}/FL/
mkdir -p ${run_path}/CCS/

# Full length spliting
# blastn: locating primers

/software/python/envs/py38/bin/samtools view $bam | awk '{print ">"$1"\n"$10}' > ${run_path}/CCS/ccs.fa
/software/ncbi-blast-2.10.0+/bin/blastn -query ${run_path}/CCS/ccs.fa  -db ${primer_db} -outfmt 7 -word_size 5 -num_threads 8 > ${run_path}/CCS/ccs.m7

# umilen 12
perl ${script_path}/classify_by_primer.pl -blastm7  ${run_path}/CCS/ccs.m7  -ccsfq ${run_path}/CCS/ccs.fq  -min_primerlen 16 -min_seqlen 50 -outdir ${run_path}/FL/

