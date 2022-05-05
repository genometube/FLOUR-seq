import pandas as pd
import sys

#################### USAGE ####################
#ARGV[1] file format: ID \s+ polyT_end_position; ID \s+ polyA_start_position
#ARGV[2] file format: input fq file (with polyT/A) REQUIRE FULL PATH
#ARGV[3] file format: output fq file name (without polyT/A) REQUIRE FULL PATH
#ARGV[4] polyT/polyA
#eg. python3 cut_pT_pA.py forward_pT_clean.csv forward.fq forward_cut5.fq polyT
###############################################
flag = str(sys.argv[4])

ID_pos = pd.read_csv(str(sys.argv[1]),header=None,sep="\s+")
dict_pos = {}
for i in range (len(ID_pos)):
    ID = ID_pos.iloc[i][0]
    pos = ID_pos.iloc[i][1]
    dict_pos[ID] = pos

fq_in = open(str(sys.argv[2]),"r")
fq_out = open(str(sys.argv[3]),"w")
if flag == "polyT":
    noT = open("noT.fq","w")
if flag == "polyA":
    noA = open("noA.fq","w")
while True:
    name = fq_in.readline().strip()
    ID = name.partition(" runid")[0]
    if not name:
        break
    seq = fq_in.readline().strip()
    strand = fq_in.readline().strip()
    qual = fq_in.readline().strip()
    
    if flag == "polyT":
        if ID in dict_pos:
            fq_out.write(name+"\n")
            fq_out.write(seq[dict_pos[ID]:]+"\n")
            fq_out.write(strand+"\n")
            fq_out.write(qual[dict_pos[ID]:]+"\n")
        else:
            noT.write(name+"\n")
            noT.write(seq+"\n")
            noT.write(strand+"\n")
            noT.write(qual+"\n")

    if flag == "polyA":
        if ID in dict_pos:
            fq_out.write(name+"\n")
            fq_out.write(seq[:dict_pos[ID]]+"\n")
            fq_out.write(strand+"\n")
            fq_out.write(qual[:dict_pos[ID]]+"\n")
        else:
            noA.write(name+"\n")
            noA.write(seq+"\n")
            noA.write(strand+"\n")
            noA.write(qual+"\n")
fq_in.close()
fq_out.close()
if flag == "polyT":
    noT.close()
if flag == "polyA":
    noA.close()
