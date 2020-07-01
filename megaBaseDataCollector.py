import sys, getopt, gzip, csv, fileinput

class megabaseInfo:
    def __init__(self, start, end, count):
        self.start = start
        self.end = end
        self.count = count


def megabaseCountMerge(file, overlap, outputPrefix):
    megabaseSize = 1000000
    chromInfoDict = {}
    currentChrom = ""
    currentMegaBaseIndex = 0
    for line in fileinput.input('-'):
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
            chromInfoDict[s[0]] = [
                megabaseInfo(0, 0 + megabaseSize, 0)]
            currentMegaBaseIndex = 0
            currentMegaBaseEnd = 0 + megabaseSize
            currentChrom = s[0]

        while (int(s[1]) > currentMegaBaseEnd):
            newMegaBaseStart = currentMegaBaseEnd - \
                (overlap * megabaseSize)
            chromInfoDict[s[0]].append(megabaseInfo(
                newMegaBaseStart, newMegaBaseStart + megabaseSize, 0))

            currentMegaBaseEnd = newMegaBaseStart + megabaseSize
            currentMegaBaseIndex += 1

        if (int(s[1]) <= currentMegaBaseEnd and int(s[1]) > currentMegaBaseEnd - (megabaseSize * overlap)):
            if currentMegaBaseIndex == len(chromInfoDict[s[0]]) - 1:
                newMegaBaseStart = currentMegaBaseEnd - \
                    (overlap * megabaseSize)
                chromInfoDict[s[0]].append(megabaseInfo(
                    newMegaBaseStart, newMegaBaseStart + megabaseSize, 1))
                chromInfoDict[s[0]][currentMegaBaseIndex].count += 1
            else:
                chromInfoDict[s[0]][currentMegaBaseIndex].count += 1
                chromInfoDict[s[0]
                                ][currentMegaBaseIndex + 1].count += 1
        else:
            chromInfoDict[s[0]][currentMegaBaseIndex].count += 1

    wCounts = csv.writer(open(outputPrefix + "megaBaseCounts.csv", "w"))
    for key in chromInfoDict.keys():
        for val in chromInfoDict[key]:
            wCounts.writerow([key, val])

    return chromInfoDict

def megabaseCount(file, overlap, outputPrefix):
    megabaseSize = 1000000
    chromInfoDict = {}
    currentChrom = ""
    currentMegaBaseIndex = 0
    if(file[-3:] == ".gz"):
        with gzip.open(file, mode='rt') as f:
            for line in f:
                if len(line.strip()) == 0:
                    continue
                elif (line[:2] == "##"):
                    continue
                elif (line[:1] == "#" ):
                    continue
                s = line.strip().split('\t')
                #New chromosome means we add a new key to the dictionary and append a new megabase counter.
                ##s[0] is the chromosome, s[1] is the position
                if (s[0] != currentChrom):
                    chromInfoDict[s[0]] = [megabaseInfo(0, 0 + megabaseSize, 0)]
                    currentMegaBaseIndex = 0
                    currentMegaBaseEnd = 0 + megabaseSize
                    currentChrom = s[0]

                while (s[1] > currentMegaBaseEnd):
                    newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
                    chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, 0))

                    currentMegaBaseEnd = newMegaBaseStart + megabaseSize
                    currentMegaBaseIndex += 1
                    
                if (s[1] <= currentMegaBaseEnd and s[1] > currentMegaBaseEnd - (megabaseSize * overlap)):
                    if currentMegaBaseIndex == len(chromInfoDict[s[0]]) -1:
                        newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
                        chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, 1))
                        chromInfoDict[s[0]][currentMegaBaseIndex].count += 1
                    else:
                        chromInfoDict[s[0]][currentMegaBaseIndex].count += 1
                        chromInfoDict[s[0]][currentMegaBaseIndex + 1].count += 1
                else:
                    chromInfoDict[s[0]][currentMegaBaseIndex].count += 1

    wCounts = csv.writer(open(outputPrefix + "megaBaseCounts.csv", "w"))
    for key in chromInfoDict.keys():
        for val in chromInfoDict[key]:
            wCounts.writerow([key,val])

    return chromInfoDict


def main(argv):
    opts, args = getopt.getopt(argv, "hp:", ['merge=', 'output=', 'overlap='])
    outputPrefix = ""
    path = ""

    for opt, arg in opts:
        if opt == '--merge':
            if arg == 'True':
                merge = True
            else:
                merge = False
        elif opt == '-p':
            path = arg
        elif opt == '--output':
            outputPrefix = arg
        elif opt == '--overlap':
            overlap = float(arg)
        elif opt in ('-h'):
            print("Use", ['--merge', '--output', '--overlap'], "to adjust parameters")
            sys.exit(0)

    
    megabaseCountMerge(path, overlap, outputPrefix)
    

if __name__ == '__main__':
      main(sys.argv[1:])
