#! /ifswh1/BC_PUB/biosoft/pipeline/10x/RNA/10x_RNA_2017a/cellranger-3.1.0/miniconda-cr-cs/4.3.21-miniconda-cr-cs-c10/bin/python

import argparse
parser = argparse.ArgumentParser(description='Extract barcode quality, generate temporary barcode table and in-whitelist barcode prior distribution.')
parser.add_argument('-i','--input', required=True, help = 'Barcodes of FL reads in fastq format.')
parser.add_argument('-w', '--whitelist', required=True, help = 'Path of the barcode whitelist.')
parser.add_argument('-t', '--threshold', required=True,type = float,help = 'Threshold of p-value for testing whether a barcode should be corrected.')
parser.add_argument('--cellranger_path', default = '/ifswh1/BC_PUB/biosoft/pipeline/10x/RNA/10x_RNA_2017a/cellranger-3.1.0/', help = 'Main path of Cellrager.')
parser.add_argument('-o','--output', required=True, help = 'Output fasta with barcodes corrected using whitelist.')

args = parser.parse_args()

import sys
sys.path.append(args.cellranger_path+"/cellranger-cs/3.1.0/tenkit/lib/python/")
sys.path.append(args.cellranger_path+"/cellranger-cs/3.1.0/lib/python/")
import string
import array
import numpy as np
import cellranger.stats as cr_stats
import tenkit.constants as tk_constants
import tenkit.seq as tk_seq
import gzip
from collections import Counter

if args.whitelist[-3:] == '.gz':
    wl = set(gzip.open(args.whitelist).read().split('\n'))
else:
    wl = set(open(args.whitelist).read().split('\n'))
wl_dict = Counter()
bc_dict_nwl = {}
my_counter = [0,0,0]

with open(args.input,'r') as f,open(args.output+'barcode_valid.fa','w') as g:
    line = f.readline()
    while line:
        readid = line.strip()[1:]
        barcode = tk_seq.get_rev_comp(f.readline().strip())
        f.readline()
        qual = f.readline().strip()[::-1]
        if barcode in wl:
            wl_dict.update([barcode])
            g.write(''.join(['>'+readid,'\n',barcode,'\n']))
            my_counter[0] += 1
        else:
            bc_dict_nwl[readid] = [barcode,qual]
            my_counter[1] += 1
        line = f.readline()

wl_sum = np.array(list(wl_dict.values())).sum()
wl_dict = dict(wl_dict)
for i in wl_dict:
    wl_dict[i] = float(wl_dict[i])/wl_sum

with open(args.output+'barcode_valid.fa','a') as f:
    for i in bc_dict_nwl:
        crbc = cr_stats.correct_bc_error(args.threshold,bc_dict_nwl[i][0],bc_dict_nwl[i][1],wl_dict)
        if crbc:
            f.write(''.join(['>'+i,'\n',crbc,'\n']))
            my_counter[2] += 1

with open(args.output+'barcode_correction.log','w') as f:
    f.write('\t'.join(['Barcode','Barcode_in_whitelist','Barcode_notin_whitelist','Barcode_corrected','Barcode_valid'])+'\n')
    f.write('\t'.join([str(my_counter[0]+my_counter[1]),str(my_counter[0]),str(my_counter[1]),str(my_counter[2]),str(my_counter[0]+my_counter[2])])+'\n')
