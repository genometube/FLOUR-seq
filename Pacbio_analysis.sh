script_path=/script

### call CCS
sh ${script_path}/Call_CCS.sh

### Classify CCS 
sh ${script_path}/Classify_CCS.sh

### Barcode correction 
sh  ${script_path}/barcode_correction.sh

### Mapping
sh ${script_path}/mapping.sh


### Mapping reads to reference annotation
sh ${script_path}/gffcompare.sh

### Extract exon intron
sh ${script_path}/exon_intron.sh

### Generate gene-barcode exon intron-cell matrix
sh ${script_path}/expr.sh

### cell-calling using gene-barcode matrix
sh ${script_path}/cell-calling.sh
