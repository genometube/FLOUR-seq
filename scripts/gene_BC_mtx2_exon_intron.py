#*coding=utf-8
import optparse,os,sys
import re
import argparse
import collections
from operator import itemgetter

usage='''
Descript:
Author:chenweitian
Email:chenweitian@bgi.com
Date:20210803
'''
option = optparse.OptionParser(usage)
option.add_option('','--input',help='Tmap file generated from gffcompare',default='' )
option.add_option('','--output',help='',default='Path to output matirx.' )
option.add_option('','--barcode',help='',default='Corrected barcodes of FL reads in fastq format.' )
option.add_option('','--umi',help='',default='Corrected UMI of FL reads in fasta format.' )
(opts, args) = option.parse_args()

umi_opt = opts.umi
barcode_opt = opts.barcode
output_opt = opts.output
input_opt = opts.input


def main():


    ###############################################3
    id2barcode_umi = {}
    print ("gene\tbarcode\tcount\texon\tintron\tu\ts")
    with open(umi_opt,'r') as f:
        #d = [f.readline().strip() for i in range(4)]
        d = [f.readline().strip() for i in range(2)]
        while d[0]:
            id2barcode_umi[d[0][1:]] = ['',d[1]]
            #d = [f.readline().strip() for i in range(4)]
            d = [f.readline().strip() for i in range(2)]
    
    with open(barcode_opt,'r') as f:
        #d = [f.readline().strip() for i in range(2)]
        d = [f.readline().strip() for i in range(4)]
        while d[0]:
            if id2barcode_umi.get(d[0][1:],0) != 0:
                id2barcode_umi[d[0][1:]][0] = d[1]
            #d = [f.readline().strip() for i in range(2)]
            d = [f.readline().strip() for i in range(4)]

    #4b2f-4ae3-89e2-95fac6ab87e0       ['cell_label_885082', 'umi_885082']
    bu2gene = {}
    exon_intron = {}
    deumi = {}
    with open(input_opt,'r') as f:
        for line in f:
            if line.startswith("ref_gene_id"):continue
            d = line.strip().split('\t')
            if (d[0] != '-') and (id2barcode_umi.get(d[4],0) != 0):
                key = (id2barcode_umi[d[4]][0],id2barcode_umi[d[4]][1],d[0])
                key2 = (id2barcode_umi[d[4]][0],d[0])

                if  key  not in deumi:
                    deumi[key] = 1
                else:
                    continue
                if key2 not in bu2gene:
                    #bu2gene.setdefault((id2barcode_umi[d[4]][0],id2barcode_umi[d[4]][1]),{})[d[1]]=1
                    bu2gene[key2] = 1
                    exon_intron.setdefault(key2,{})['exon'] = int(d[-2])
                    exon_intron.setdefault(key2,{})['intron'] = int(d[-1])
                    if int(d[-1]) > 0:
                        exon_intron.setdefault(key2,{})['u'] = 1
                        exon_intron.setdefault(key2,{})['s'] = 0
                    else:
                        exon_intron.setdefault(key2,{})['u'] = 0
                        exon_intron.setdefault(key2,{})['s'] = 1

                else:
                    #bu2gene.setdefault((id2barcode_umi[d[4]][0],id2barcode_umi[d[4]][1]),{})[d[1]] += 1
                    bu2gene[key2] += 1
                    exon_intron.setdefault(key2,{})['exon'] += int(d[-2])
                    exon_intron.setdefault(key2,{})['intron'] += int(d[-1])
                    if int(d[-1]) > 0:
                        exon_intron.setdefault(key2,{})['u'] += 1
                        exon_intron.setdefault(key2,{})['s'] += 0
                    else:
                        exon_intron.setdefault(key2,{})['u'] += 0
                        exon_intron.setdefault(key2,{})['s'] += 1


    del id2barcode_umi
    del deumi
    for key in bu2gene:
        #b,u,gene =  key
        count = bu2gene[key]
        exon=exon_intron[key]['exon']
        intron=exon_intron[key]['intron']
        u = exon_intron[key]['u']
        s = exon_intron[key]['s']
        print ("%s\t%s\t%s\t%s\t%s\t%s\t%s"%(key[1],key[0],count,exon,intron,u,s))

if __name__ == '__main__':
    if len(sys.argv)<2:
        os.system("python %s -h"%(sys.argv[0]))
        sys.exit(1)
    else:
        main()
