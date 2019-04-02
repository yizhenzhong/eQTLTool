import sys
from bisect import bisect
def getStart(f):
        start = {}
        for line in open(f):
                items = line.split()
                try:
                        chrs =  int(items[0][3:])
                        if chrs in start:
                                start[chrs].append(int(items[1]))
                        else:
                                start[chrs] = [int(items[1])]
                except:
                        continue
        return start

def checkSorted(l):
        if l == sorted(l):
                return True
        else:
                return False

def minDist(start, chrs, pos):
        #Use positive distance for upstream SNP and negative distance for downstream SNP
        tract = start[chrs]
        index = bisect(tract, pos)
        if index == len(tract):
                return tract[-1]-pos
        if abs(tract[index-1] - pos) > abs(tract[index] - pos):
                return tract[index] -  pos
        else:
                return tract[index-1] - pos

def readBim(bim, start):
        out = open(bim[:-3]+"distTSS.txt","w")
        for line in open(bim):
                items = line.split()
                out.write("\t".join([items[0], items[1], str(minDist(start, int(items[0]), int(items[3])))])+"\n")
        out.close()                

def main():
        bed = sys.argv[1]
        start = getStart(bed)
        for chr in start:
                print(chr)
                if not checkSorted(start[chr]):
                        print("{} not sorted".format(str(chr)))
        readBim(sys.argv[2], start)


if __name__ == '__main__':
    main()
                        
