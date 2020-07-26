import sys, getopt, gzip, csv, fileinput, os

class megabaseInfo:
    def __init__(self, start, end, count, insertion, deletion, snv):
        self.start = start
        self.end = end
        self.count = count
        self.insertion = insertion
        self.deletion = deletion
        self.snv = snv

class familyInfo:
    def __init__(self, familyID, probandID, siblingID, probandGender, siblingGender, fatherAge):
        self.familyID = familyID
        self.probandID = probandID
        self.siblingID = siblingID
        self.probandGender = probandGender
        self.siblingGender = siblingGender
        self.fatherAge = fatherAge
        

def pairSiblings(famFile, sampleFile):
    families = []

    with open(sampleFile, mode='rt') as f:
        header = f.readline()
        for line in f:
            s = line.strip().split("\t")
            families.append(familyInfo(s[0], s[3], s[4], "?", "?", s[19]))

    familyIndex = 0
    currentFamilyId = 0

    with open(famFile, mode='rt') as f:
        header = f.readline()
        for line in f:
            s = line.strip().split()
            ind = s[1].split(".")[1]
            if int(s[0]) != int(families[familyIndex].familyID):
                while int(s[0]) > int(families[familyIndex].familyID):
                    #start new family set.
                    familyIndex += 1
                if (int(s[0]) < int(families[familyIndex].familyID)):
                    continue
            else:
                if ind == "p1":
                    families[familyIndex].probandGender = s[3]
                elif ind == "s1":
                    families[familyIndex].siblingGender = s[3]
                
                if families[familyIndex].siblingGender != "?" and families[familyIndex].probandGender != "?":
                    familyIndex += 1

    return families


def appendChromInfoDict(data, megaBaseStart, megaBaseSize, chrom, variantCount, insert, Del, snv):
    data[chrom].append(megabaseInfo(megaBaseStart, megaBaseStart + megaBaseSize, variantCount, insert, Del, snv))

def updateChromInfoDict(data, index, chrom, variantCount, insert, Del, snv, updateNext):
    data[chrom][index].count += variantCount
    data[chrom][index].insertion += insert
    data[chrom][index].deletion += Del
    data[chrom][index].snv += snv
    if updateNext:
        data[chrom][index + 1].count += variantCount
        data[chrom][index + 1].insertion += insert
        data[chrom][index + 1].deletion += Del
        data[chrom][index + 1].snv += snv

            
def megabaseCountMergeFamily(file, overlap, binsize, outputPrefix, familyData):
    megabaseSize = binsize
    probandDataSet = True
    chromInfoDict = {}
    chromInfoDictMM = {}
    chromInfoDictMF = {}
    chromInfoDictFM = {}
    chromInfoDictFF = {}
    #[MM, MF, FM, FF]
    result = [{}, {}, {}, {}]

    for root, dirs, files in os.walk(file):
        temproot = root.strip().split("/")
        print("root", temproot)
        if temproot[len(temproot) - 1] == "s1":
            probandDataSet = False
        for filename in files:
            if(filename[-13:] == ".FINAL.vcf.gz"):
                with gzip.open(root + filename, mode='rt') as f:
                    print("Scanning variants from file:", filename)
                    currentChrom = ""
                    currentMegaBaseIndex = 0
                    currentDataSet = 0
                    for line in f:
                        if len(line.strip()) == 0:
                            continue
                        elif (line[:2] == "##"):
                            continue
                        elif (line[:1] == "#"):
                            #parse which sample the current file has.
                            sample = line.strip().split('\t')[9][2:]
                            print("Sample:", sample)
                            if probandDataSet:
                                for x in familyData:
                                    if x.probandID == sample:
                                        if (x.probandGender == "male") and (x.siblingGender == "male"):
                                            currentDataSet = 0
                                        elif (x.probandGender == "male") and (x.siblingGender == "female"):
                                            currentDataSet = 1
                                        elif (x.probandGender == "female") and (x.siblingGender == "male"):
                                            currentDataSet = 2
                                        elif (x.probandGender == "female") and (x.siblingGender == "female"):
                                            currentDataSet = 3
                                        break            
                            else:
                                for x in familyData:
                                    if x.siblingID == sample:
                                        if (x.probandGender == "male") and (x.siblingGender == "male"):
                                            currentDataSet = 0
                                        elif (x.probandGender == "male") and (x.siblingGender == "female"):
                                            currentDataSet = 1
                                        elif (x.probandGender == "female") and (x.siblingGender == "male"):
                                            currentDataSet = 2
                                        elif (x.probandGender == "female") and (x.siblingGender == "female"):
                                            currentDataSet = 3
                                        break    
                            continue

                        s = line.strip().split('\t')
                        IDField = s[2].split('-')[0]

                        insert = 0
                        Del = 0
                        snv = 0

                        #Parse if the variant is an insertion, deletion or both
                        if any(x in IDField for x in ["BND", '+', '_', "DUP"]):
                            continue
                        else:
                            if ("I" in IDField):
                                insert += 1
                            if ("Y" in IDField):
                                insert += 1
                            if ("D" in IDField):
                                Del += 1
                            if ("X" in IDField):
                                snv += 1

                        variantCount = insert + Del + snv

                        #Place the variant in the right megabase bin
                        #New chromosome means we add a new key to the dictionary and append a new megabase counter.
                        ##s[0] is the chromosome, s[1] is the position
                        if (len(s[0]) > 5):
                            continue
                        if (s[0] != currentChrom):
                            # if (s[0] not in chromInfoDict):
                            #     chromInfoDict[s[0]] = [megabaseInfo(0, 0 + megabaseSize, 0, 0, 0, 0)]    
                            #     currentMegaBaseEnd = 0 + megabaseSize
                            if (s[0] not in result[currentDataSet]):
                                result[currentDataSet][s[0]] = [megabaseInfo(0, 0 + megabaseSize, 0, 0, 0, 0)]    
                                currentMegaBaseEnd = 0 + megabaseSize
                            
                            currentMegaBaseIndex = 0
                            currentMegaBaseEnd = result[currentDataSet][s[0]][0].end
                            currentChrom = s[0]

                        while (int(s[1]) > currentMegaBaseEnd):
                            #If the megabase isnt initialized yet
                            if (currentMegaBaseIndex + 1 == len(result[currentDataSet][s[0]])):
                                newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
                                #chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, 0, 0, 0, 0))
                                appendChromInfoDict(result[currentDataSet], newMegaBaseStart, megabaseSize, s[0], 0, 0, 0, 0)
                                currentMegaBaseEnd = newMegaBaseStart + megabaseSize
                                currentMegaBaseIndex += 1
                            else:
                                currentMegaBaseIndex += 1
                                currentMegaBaseEnd = result[currentDataSet][s[0]][currentMegaBaseIndex].end
                                

                        #If we find a variant that falls in the overlap for two megabases.
                        if (int(s[1]) <= currentMegaBaseEnd and int(s[1]) > currentMegaBaseEnd - (megabaseSize * overlap)):
                            if (currentMegaBaseIndex == len(result[currentDataSet][s[0]]) - 1):
                                newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
                                appendChromInfoDict(result[currentDataSet], newMegaBaseStart, megabaseSize, s[0], variantCount, insert, Del, snv)
                                #chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, variantCount, insert, Del, snv))
                                updateChromInfoDict(result[currentDataSet], currentMegaBaseIndex, s[0], variantCount, insert, Del, snv, False)
                                # chromInfoDict[s[0]][currentMegaBaseIndex].count += variantCount
                                # chromInfoDict[s[0]][currentMegaBaseIndex].insertion += insert
                                # chromInfoDict[s[0]][currentMegaBaseIndex].deletion += Del
                                # chromInfoDict[s[0]][currentMegaBaseIndex].snv += snv
                            else:
                                updateChromInfoDict(result[currentDataSet], currentMegaBaseIndex, s[0], variantCount, insert, Del, snv, True)
                                # chromInfoDict[s[0]][currentMegaBaseIndex].count += variantCount
                                # chromInfoDict[s[0]][currentMegaBaseIndex].insertion += insert
                                # chromInfoDict[s[0]][currentMegaBaseIndex].deletion += Del
                                # chromInfoDict[s[0]][currentMegaBaseIndex].snv += snv
                                # chromInfoDict[s[0]][currentMegaBaseIndex + 1].count += variantCount
                                # chromInfoDict[s[0]][currentMegaBaseIndex + 1].insertion += insert
                                # chromInfoDict[s[0]][currentMegaBaseIndex + 1].deletion += Del
                                # chromInfoDict[s[0]][currentMegaBaseIndex + 1].snv += snv
                        else:
                            updateChromInfoDict(result[currentDataSet], currentMegaBaseIndex, s[0], variantCount, insert, Del, snv, False)
                            # chromInfoDict[s[0]][currentMegaBaseIndex].count += variantCount
                            # chromInfoDict[s[0]][currentMegaBaseIndex].insertion += insert
                            # chromInfoDict[s[0]][currentMegaBaseIndex].deletion += Del
                            # chromInfoDict[s[0]][currentMegaBaseIndex].snv += snv
                            
                        #print(chrom, int(s[1]), chromInfoDict[chrom][currentMegaBaseIndex].end - megabaseSize, chromInfoDict[chrom][currentMegaBaseIndex].end, chromInfoDict[chrom][currentMegaBaseIndex].count)

    outputIndex = 0
    typePrefix = ["MM", "MF", "FM", "FF"]
    for x in result:
        wCounts = csv.writer(open(outputPrefix + "." + typePrefix[outputIndex] + ".megaBaseCounts.csv", "w"))
        wCounts.writerow(["Chrom", "Start", "End", "Count", "Insertions", "Deletions", "SNV"])
        for key in x.keys():
            for val in x[key]:
                wCounts.writerow([key, val.start, val.end, val.count, val.insertion, val.deletion, val.snv])

    return result




def megabaseCountMerge(file, overlap, binsize, outputPrefix):
    megabaseSize = binsize
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
                        IDField = s[2].split('-')[0]

                        insert = 0
                        Del = 0
                        snv = 0

                        #Parse if the variant is an insertion, deletion or both
                        if any(x in IDField for x in ["BND", '+', '_', "DUP"]):
                            continue
                        else:
                            if ("I" in IDField):
                                insert += 1
                            if ("Y" in IDField):
                                insert += 1
                            if ("D" in IDField):
                                Del += 1
                            if ("X" in IDField):
                                snv += 1

                        variantCount = insert + Del + snv

                        #Place the variant in the right megabase bin
                        #New chromosome means we add a new key to the dictionary and append a new megabase counter.
                        ##s[0] is the chromosome, s[1] is the position
                        if (len(s[0]) > 5):
                            continue
                        if (s[0] != currentChrom):
                            if (s[0] not in chromInfoDict):
                                chromInfoDict[s[0]] = [megabaseInfo(0, 0 + megabaseSize, 0, 0, 0, 0)]    
                                currentMegaBaseEnd = 0 + megabaseSize
                            
                            currentMegaBaseIndex = 0
                            currentMegaBaseEnd = chromInfoDict[s[0]][0].end
                            currentChrom = s[0]

                        while (int(s[1]) > currentMegaBaseEnd):
                            #If the megabase isnt initialized yet
                            if (currentMegaBaseIndex + 1 == len(chromInfoDict[s[0]])):
                                newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
                                chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, 0, 0, 0, 0))
                                currentMegaBaseEnd = newMegaBaseStart + megabaseSize
                                currentMegaBaseIndex += 1
                            else:
                                currentMegaBaseIndex += 1
                                currentMegaBaseEnd = chromInfoDict[s[0]][currentMegaBaseIndex].end
                                

                        #If we find a variant that falls in the overlap for two megabases.
                        if (int(s[1]) <= currentMegaBaseEnd and int(s[1]) > currentMegaBaseEnd - (megabaseSize * overlap)):
                            if (currentMegaBaseIndex == len(chromInfoDict[s[0]]) - 1):
                                newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
                                chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, variantCount, insert, Del, snv))
                                chromInfoDict[s[0]][currentMegaBaseIndex].count += variantCount
                                chromInfoDict[s[0]][currentMegaBaseIndex].insertion += insert
                                chromInfoDict[s[0]][currentMegaBaseIndex].deletion += Del
                                chromInfoDict[s[0]][currentMegaBaseIndex].snv += snv
                            else:
                                chromInfoDict[s[0]][currentMegaBaseIndex].count += variantCount
                                chromInfoDict[s[0]][currentMegaBaseIndex].insertion += insert
                                chromInfoDict[s[0]][currentMegaBaseIndex].deletion += Del
                                chromInfoDict[s[0]][currentMegaBaseIndex].snv += snv
                                chromInfoDict[s[0]][currentMegaBaseIndex + 1].count += variantCount
                                chromInfoDict[s[0]][currentMegaBaseIndex + 1].insertion += insert
                                chromInfoDict[s[0]][currentMegaBaseIndex + 1].deletion += Del
                                chromInfoDict[s[0]][currentMegaBaseIndex + 1].snv += snv
                        else:
                            chromInfoDict[s[0]][currentMegaBaseIndex].count += variantCount
                            chromInfoDict[s[0]][currentMegaBaseIndex].insertion += insert
                            chromInfoDict[s[0]][currentMegaBaseIndex].deletion += Del
                            chromInfoDict[s[0]][currentMegaBaseIndex].snv += snv
                            
                        #print(chrom, int(s[1]), chromInfoDict[chrom][currentMegaBaseIndex].end - megabaseSize, chromInfoDict[chrom][currentMegaBaseIndex].end, chromInfoDict[chrom][currentMegaBaseIndex].count)

    wCounts = csv.writer(open(outputPrefix + "megaBaseCounts.csv", "w"))
    wCounts.writerow(["Chrom", "Start", "End", "Count", "Insertions", "Deletions", "SNV"])
    for key in chromInfoDict.keys():
        for val in chromInfoDict[key]:
            wCounts.writerow([key, val.start, val.end, val.count, val.insertion, val.deletion, val.snv])

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
        chrom = s[0].split('_')[0]
        IDField = s[2].split('-')[0]

        insert = 0
        Del = 0
        snv = 0

        #Parse if the variant is an insertion, deletion or both
        if any(x in IDField for x in ["BND", '+', '_', "DUP"]):
            continue
        else:
            if ("I" in IDField):
                insert += 1
            if ("Y" in IDField):
                insert += 1
            if ("D" in IDField):
                Del += 1
            if ("X" in IDField):
                snv += 1

        variantCount = insert + Del + snv

        #Place the variant in the right megabase bin
        if (chrom != currentChrom):
            print("Scanning variants in:", chrom)
            if (chrom not in chromInfoDict):
                chromInfoDict[s[0]] = [megabaseInfo(0, 0 + megabaseSize, 0, 0, 0, 0)]    
                currentMegaBaseEnd = 0 + megabaseSize
            
            currentMegaBaseEnd = chromInfoDict[s[0]][0].end
            currentMegaBaseIndex = 0
            currentChrom = chrom
            

        while (int(s[1]) > currentMegaBaseEnd):
            #If the megabase isnt initialized yet
            if (currentMegaBaseIndex + 1 == len(chromInfoDict[s[0]]) ):
                newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
                chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, 0, 0, 0 ,0))
                currentMegaBaseEnd = newMegaBaseStart + megabaseSize
                currentMegaBaseIndex += 1
            else:
                currentMegaBaseIndex += 1
                currentMegaBaseEnd = chromInfoDict[s[0]][currentMegaBaseIndex].end
                

        #If we find a variant that falls in the overlap for two megabases.
        if (int(s[1]) <= currentMegaBaseEnd and int(s[1]) > currentMegaBaseEnd - (megabaseSize * overlap)):
            if (currentMegaBaseIndex == len(chromInfoDict[s[0]]) - 1):
                newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
                chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, variantCount, insert, Del, snv))
                chromInfoDict[s[0]][currentMegaBaseIndex].count += variantCount
                chromInfoDict[s[0]][currentMegaBaseIndex].insertion += insert
                chromInfoDict[s[0]][currentMegaBaseIndex].deletion += Del
                chromInfoDict[s[0]][currentMegaBaseIndex].snv += snv
            else:
                chromInfoDict[s[0]][currentMegaBaseIndex].count += variantCount
                chromInfoDict[s[0]][currentMegaBaseIndex].insertion += insert
                chromInfoDict[s[0]][currentMegaBaseIndex].deletion += Del
                chromInfoDict[s[0]][currentMegaBaseIndex].snv += snv
                chromInfoDict[s[0]][currentMegaBaseIndex + 1].count += variantCount
                chromInfoDict[s[0]][currentMegaBaseIndex + 1].insertion += insert
                chromInfoDict[s[0]][currentMegaBaseIndex + 1].deletion += Del
                chromInfoDict[s[0]][currentMegaBaseIndex + 1].snv += snv
        else:
            chromInfoDict[s[0]][currentMegaBaseIndex].count += variantCount
            chromInfoDict[s[0]][currentMegaBaseIndex].insertion += insert
            chromInfoDict[s[0]][currentMegaBaseIndex].deletion += Del
            chromInfoDict[s[0]][currentMegaBaseIndex].snv += snv

    wCounts = csv.writer(open(outputPrefix + "megaBaseCounts.csv", "w"))
    wCounts.writerow(["Chrom", "Start", "End", "Count", "Insertions", "Deletions", "SNV"])
    for key in chromInfoDict.keys():
        for val in chromInfoDict[key]:
            wCounts.writerow([key, val.start, val.end, val.count, val.insertion, val.deletion])

    return chromInfoDict


def main(argv):
    opts, args = getopt.getopt(argv, "h p:", ['merge=', 'output=', 'overlap=', 'binsize=', 'famFile=', 'sampleFile='])
    outputPrefix = ""
    path = '-'
    merge = False
    binsize = 1000000
    famFile = ""

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
        elif opt == '--famFile':
            famFile = arg
        elif opt == '--sampleFile':
            sampleFile = arg
        elif opt in ('-h'):
            print("Use", ['--merge', '--output', '--overlap'], "to adjust parameters")
            sys.exit(0)

    familyData = pairSiblings(famFile, sampleFile)

    if (merge):
        if famFile == "":
            megabaseCountMerge(path, overlap, binsize, outputPrefix)
        else:
            megabaseCountMergeFamily(path, overlap, binsize, outputPrefix, familyData)
    else:
        megabaseCount(path, overlap, binsize, outputPrefix)
    

if __name__ == '__main__':
      main(sys.argv[1:])
