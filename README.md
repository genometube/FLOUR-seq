# HISOFA-seq
Pipeline for HIgh-throughput Single-cell ONT Full-length RNA sequencing (HISOFA-seq)
## Barcode identification
The code for barcode demultiplexing step is recorded in "debarcode_HISOFA.sh". "cell_label.zip" is a zipped fasta file of barcode reference. "Find_polyA.py" and "cut_pT_pA.py" are python scripts used in the pipeline. Other tools required in this step include fastp, seqtk, bowtie2 and samtools. 
## Matrix to velocity
The script "parepare_matrix_to_velocity.sh" is used to prepare matrices for velocity analysis in R. RNA velocity can refer https://github.com/velocyto-team/velocyto.R. Region velocity can refer https://github.com/Dekayzc/Regionvelocity. 
