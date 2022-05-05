#*coding=utf-8
import optparse,os,sys
import gzip
usage='''
Descript:
Author:chenweitian
Email:chenweitian@bgi.com
Date:20210107
\n
input :fastq 
Output:read_id\tstart,end,length,poly_length
\n
Run_demo:python3 %s  --fastq urFastqFile --polyA yes 
'''%(sys.argv[0])
option = optparse.OptionParser(usage)
option.add_option('','--fastq',help='',default='' )
option.add_option('','--polyA',help='default:yes',default='yes' )
option.add_option('','--polyT',help='default:no',default='no' )

(opts, args) = option.parse_args()

fastq_opt = opts.fastq
polyA_opt = opts.polyA
polyT_opt = opts.polyT


def InvervalMerger(Intervals,MinSeed, MinLength):
    i = 0
    PolyIntervals = []
    while i < len(Intervals):
        if Intervals[i][2]>MinSeed:
            left = i - 1
            right = i + 1
            FrameLength = Intervals[i][2]
            FrameCount = Intervals[i][2]
            TmpIndex = [Intervals[i][0],Intervals[i][1]]
            flag = ((FrameCount/FrameLength<0.7) or (left==0) or (right == len(Intervals)))
            while not flag:             
                if (Intervals[left][2]/(Intervals[left+1][0]-Intervals[left][0])<0.3) and (Intervals[right][2]/(Intervals[right][1]-Intervals[right-1][1]))<0.3:
                    break
                else:
                    leftLength = FrameLength+Intervals[left+1][0]-Intervals[left][0]
                    leftCount = FrameCount+Intervals[left][2]
                    rightLength = FrameLength+Intervals[right][1]-Intervals[right-1][1]
                    rightCount = FrameCount+Intervals[right][2]
                    if leftCount/leftLength>rightCount/rightLength:
                        FrameLength = leftLength
                        FrameCount = leftCount
                        TmpIndex[0] = Intervals[left][0]
                        left -= 1
                    else:
                        FrameLength = rightLength
                        FrameCount = rightCount
                        TmpIndex[1] = Intervals[right][1]
                        right += 1
                flag = ((FrameCount/FrameLength<0.7) or (left==0) or (right == len(Intervals)))
            if FrameLength >=MinLength:
                PolyIntervals.append(TmpIndex+[FrameLength,FrameCount])
            i = right
        else:
            i += 1
    return(PolyIntervals)
    
def IntervalFinder(seq):
    leftA = 0
    flagA = 0
    leftT = 0
    flagT = 0
    intervalsA = []
    intervalsT = []
    for i in range(len(seq)):
        if flagA == 0:
            if seq[i] == 'A':
                leftA = i
                flagA = 1
        else:
            if seq[i] != 'A':
                intervalsA.append([leftA,i,i-leftA])
                flagA = 0
        if flagT == 0:
            if seq[i] == 'T':
                leftT = i
                flagT = 1
        else:
            if seq[i] != 'T':
                intervalsT.append([leftT,i,i-leftT])
                flagT = 0
    if flagA ==1:
        intervalsA.append([leftA,i,i-leftA])
    elif flagT == 1:
        intervalsT.append([leftT,i,i-leftT])
    return([intervalsA,intervalsT])
def print_fq(name,inter,seq):
    tmp = []
    for i in inter:
        #i[1]=i[1]+1
        i[0]=i[0]+1
        tmp.append(','.join(map(str,i)))
    #print ("%s\t%s\n%s"%(name,';'.join(tmp),seq))
    print ("%s\t%s"%(name,';'.join(tmp)))
def main():
    if fastq_opt.endswith('.gz'):
        fq = gzip.open(fastq_opt,'r')
    else:
        fq = open(fastq_opt,'r')
    while  1:
        polyA = []
        polyT = []
        name = fq.readline().strip()
        if not name:
            break
        seq = fq.readline().strip()
        strand = fq.readline().strip()
        qua = fq.readline().strip()        
        intervals = IntervalFinder(seq)
        if polyA_opt == "yes": 
            polyA = InvervalMerger(intervals[0],3,8)
            if len(polyA) == 0:continue
            #print polyA
            #print (seq)
            print_fq(name,polyA,seq)
        if polyT_opt == "yes":    
            polyT = InvervalMerger(intervals[1],3,8)
            if len(polyT) == 0:continue

            print_fq(name,polyT,seq)
        '''
        for site in polyA:
            print ("%s\t%s"%(name,site))
        '''
def main_test():
    seq="CGTTAAAGGCTAGGCTCACAACCAAAAATATAAGAGTTCGGTTCCCAGCACCCACGGCTGTCTCTCCAGCCACCAAAAAAAAAAAGAAAAAAAAAA"
    seq="AAAAAGGAAATAAAA"
    intervals = IntervalFinder(seq)
    polyA = InvervalMerger(intervals[0],1,1)
    print (polyA)

if __name__ == '__main__':
    if len(sys.argv)<2:
        os.system("python %s -h"%(sys.argv[0]))
        sys.exit(1)
    else:
        main()

