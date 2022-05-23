#*coding=utf-8
import sys
dicta = {}
for line in open(sys.argv[1],'r'):
    line = line.strip().split("\t")
    type = line[3]
    isoform = line[7]
    length = int(line[8])
    if length >20:
        dicta.setdefault(isoform,[]).append(line[3])

for read in dicta:
    exon = dicta[read].count('exon_DIY')
    intron = dicta[read].count('intron_DIY')
    print ("%s\t%s\t%s"%(read,exon,intron))
