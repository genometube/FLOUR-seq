# Call CCS
export PATH=$PATH:/smrtlink/smrtcmds/bin
export PATH=$PATH:/ncbi-blast-2.10.0+/bin
export PATH=$PATH:/R-3.4.1/bin

ccs *.subreads.bam ccs.bam --min-passes 0 --min-length 50 --max-length 21000 --min-rq 0.75 -j 20