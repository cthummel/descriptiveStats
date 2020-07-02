import sys, getopt, gzip, csv, fileinput, os

class megabaseInfo:
    def __init__(self, start, end, count):
        self.start = start
        self.end = end
        self.count = count


def megabaseCountMerge(file, overlap, binsize, outputPrefix):
    megabaseSize = 1000000
    chromInfoDict = {}
    for root, dirs, files in os.walk(file):
        for filename in files:
            if(filename[-13:] == ".FINAL.vcf.gz"):
                with gzip.open(root + filename, mode='rt') as f:
                    print("Scanning variants from file:", filename)
                    currentChrom = ""
                    currentMegaBaseIndex = 0
                    for line in f:
                        if len(line.strip()) == 0:
                            continue
                        elif (line[:2] == "##"):
                            continue
                        elif (line[:1] == "#"):
                            continue
                        s = line.strip().split('\t')
                        #New chromosome means we add a new key to the dictionary and append a new megabase counter.
                        ##s[0] is the chromosome, s[1] is the position
                        if (s[0] != currentChrom):
                            #print("Scanning variants in:", s[0])
                            if (s[0] not in chromInfoDict):
                                chromInfoDict[s[0]] = [megabaseInfo(0, 0 + megabaseSize, 0)]    
                                currentMegaBaseEnd = 0 + megabaseSize
                            
                            print("oldchrom:", currentChrom, "Newchrom:", s[0])
                            currentMegaBaseIndex = 0
                            currentChrom = s[0]
                            

                        while (int(s[1]) > currentMegaBaseEnd):
                            #If the megabase isnt initialized yet
                            if (currentMegaBaseIndex + 1 == len(chromInfoDict[s[0]]) ):
                                print("new megabase")
                                newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
                                chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, 0))
                                currentMegaBaseEnd = newMegaBaseStart + megabaseSize
                                currentMegaBaseIndex += 1
                            else:
                                currentMegaBaseIndex += 1
                                currentMegaBaseEnd = chromInfoDict[s[0]][currentMegaBaseIndex].end
                                print("looking for megabase", currentMegaBaseEnd -megabaseSize, currentMegaBaseEnd, currentMegaBaseIndex)
                                

                        #If we find a variant that falls in the overlap for two megabases.
                        if (int(s[1]) <= currentMegaBaseEnd and int(s[1]) > currentMegaBaseEnd - (megabaseSize * overlap)):
                            if (currentMegaBaseIndex == len(chromInfoDict[s[0]]) - 1):
                                newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
                                chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, 1))
                                chromInfoDict[s[0]][currentMegaBaseIndex].count += 1
                            else:
                                print("count before:", chromInfoDict[s[0]][currentMegaBaseIndex].count)
                                chromInfoDict[s[0]][currentMegaBaseIndex].count += 1
                                print("count after:", chromInfoDict[s[0]][currentMegaBaseIndex].count)
                                chromInfoDict[s[0]][currentMegaBaseIndex + 1].count += 1
                        else:
                            chromInfoDict[s[0]][currentMegaBaseIndex].count += 1
                            
                        print(s[0], int(s[1]), chromInfoDict[s[0]][currentMegaBaseIndex-1].end - megabaseSize, chromInfoDict[s[0]][currentMegaBaseIndex-1].end, chromInfoDict[s[0]][currentMegaBaseIndex].count)

    wCounts = csv.writer(open(outputPrefix + "megaBaseCounts.csv", "w"))
    wCounts.writerow(["Chrom", "Start", "End", "Count"])
    for key in chromInfoDict.keys():
        for val in chromInfoDict[key]:
            wCounts.writerow([key, val.start, val.end, val.count])

    return chromInfoDict

def megabaseCount(file, overlap, binsize, outputPrefix):
    megabaseSize = 1000000
    chromInfoDict = {}
    currentChrom = ""
    currentMegaBaseIndex = 0
    for line in fileinput.input(file):
        if len(line.strip()) == 0:
            continue
        elif (line[:2] == "##"):
            continue
        elif (line[:1] == "#"):
            continue
        s = line.strip().split('\t')
        #New chromosome means we add a new key to the dictionary and append a new megabase counter.
        ##s[0] is the chromosome, s[1] is the position
        if (s[0] != currentChrom):
            print("Scanning variants in:", s[0])
            if (s[0] not in chromInfoDict):
                chromInfoDict[s[0]] = [megabaseInfo(0, 0 + megabaseSize, 0)]    
                currentMegaBaseEnd = 0 + megabaseSize
            
            currentMegaBaseIndex = 0
            currentChrom = s[0]
            

        while (int(s[1]) > currentMegaBaseEnd):
            #If the megabase isnt initialized yet
            if (currentMegaBaseIndex + 1 == len(chromInfoDict[s[0]]) ):
                newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
                chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, 0))
                currentMegaBaseEnd = newMegaBaseStart + megabaseSize
                currentMegaBaseIndex += 1
            else:
                currentMegaBaseIndex += 1
                currentMegaBaseEnd = chromInfoDict[s[0]][currentMegaBaseIndex].end
                

        #If we find a variant that falls in the overlap for two megabases.
        if (int(s[1]) <= currentMegaBaseEnd and int(s[1]) > currentMegaBaseEnd - (megabaseSize * overlap)):
            if (currentMegaBaseIndex == len(chromInfoDict[s[0]]) - 1):
                newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
                chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, 1))
                chromInfoDict[s[0]][currentMegaBaseIndex].count += 1
            else:
                chromInfoDict[s[0]][currentMegaBaseIndex].count += 1
                chromInfoDict[s[0]][currentMegaBaseIndex + 1].count += 1
        else:
            chromInfoDict[s[0]][currentMegaBaseIndex].count += 1

    wCounts = csv.writer(open(outputPrefix + "megaBaseCounts.csv", "w"))
    wCounts.writerow(["Chrom", "Start", "End", "Count"])
    for key in chromInfoDict.keys():
        for val in chromInfoDict[key]:
            wCounts.writerow([key, val.start, val.end, val.count])

    return chromInfoDict


def main(argv):
    opts, args = getopt.getopt(argv, "h p:", ['merge=', 'output=', 'overlap=', 'binsize='])
    outputPrefix = ""
    path = '-'
    merge = False
    binsize = 1000000

    for opt, arg in opts:
        if opt == '--merge':
            merge = True
        elif opt == '-p':
            path = arg
        elif opt == '--output':
            outputPrefix = arg
        elif opt == '--overlap':
            overlap = float(arg)
        elif opt == '--binsize':
            binsize = arg
        elif opt in ('-h'):
            print("Use", ['--merge', '--output', '--overlap'], "to adjust parameters")
            sys.exit(0)

    if (merge):
        megabaseCountMerge(path, overlap, binsize, outputPrefix)
    else:
        megabaseCount(path, overlap, binsize, outputPrefix)
    

if __name__ == '__main__':
      main(sys.argv[1:])
