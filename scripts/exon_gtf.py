#*coding=utf-8
#by chenweitian@bgi.com 
# 20210721
import sys
from operator import itemgetter
import re


dicta = {}
dictb = {}
id_match = open("id_match.txt",'w')
id_match.write("gene_name\ttranscript_name\tgene_id\ttranscript_id\n")
id_match_dict = {}
for line1 in open(sys.argv[1],'r'):
    line = line1.strip().split("\t")
    site = (line[3],line[4])
    attr = line[8].split("; ")
    chr = line[0]
    strand = line[6]
    gene_name = ""
    transcript_name = ""
    gene_id = ""
    transcript_id = ""
    exon_number = ""
    for i in attr:
        if re.search("gene_id",i):
            m = re.search("gene_id \"(.*)\"",i)
            gene_id = m.group(1)
        if re.search("exon_number",i):
            m = re.search("exon_number \"(.*)\"",i)
            exon_number = m.group(1)
        if re.search("gene_name",i):
            m = re.search("gene_name \"(.*)\"",i)
            gene_name = m.group(1)
        if re.search("transcript_id",i):
            m = re.search("transcript_id \"(.*)\"",i)
            transcript_id = m.group(1)
        if re.search("transcript_name",i):
            m = re.search("transcript_name \"(.*)\"",i)
            transcript_name = m.group(1)

    id_match_key = "%s\t%s\t%s\t%s\n"%(gene_id, transcript_name, gene_id, transcript_id)
    id_match_dict[id_match_key] = 1
    
    dicta.setdefault(chr,{}).setdefault(gene_id,[]).append(site)
    #print (exon_number)
    dictb.setdefault(chr,{})[gene_id]= strand
    #print ("%s\texon_DIY\t%s"%(line[0],'\t'.join(line[2:])))
    print ("%s\texon_DIY\t%s\ttranscript_id \"%s\"; gene_id \"%s\"; exon_id \"%s\""%(line[0],'\t'.join(line[2:8]), transcript_id,gene_id,exon_number))
    #print (line1.strip()) 
'''
for i in id_match_dict:
    id_match.write("%s\n"%(i))

for chr  in  dicta:
    for  gene_name in dicta[chr]:
        data = list(set(dicta[chr][gene_name]))
        data.sort(key=itemgetter(0))
        count = 0
        intron_list = []
        strand = dictb[chr][gene_name]
        for intron in data:
            #print (count)
            count += 1
            
            s = int(intron[0])
            e = int(intron[1])
            if count==1:
                S1 = int(intron[0])
                e1 = int(intron[1]) 
            if count > 1:
                
                if s > e1:
                    intron_s = e1+1
                    intron_e = s-1
                    e1 = e
                    intron_list.append((intron_s,intron_e))
                    intron_name = "intron_%s"%(count-1)
                    if intron_s>=intron_e:continue
                    gtf_attr = "%s\tintron_DIY\texon\t%s\t%s\t.\t%s\t.\ttranscript_id \"%s\"; gene_id \"%s\"; exon_id \"%s\""%(chr,intron_s,intron_e,strand,gene_name,gene_name,intron_name)
                    print (gtf_attr)
        #print ("%s\t%s\t%s"%(chr,gene_name,data))

        #print (dictb[chr][gene_name][count])
        
        #print ("%s\t%s\t%s"%(chr,gene_name,intron_list))
'''
