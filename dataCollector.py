import numpy as np
import math, sys, getopt, gzip, os, csv, binFinder, re

class megabaseInfo:
    def __init__(self, start, end, count):
        self.start = start
        self.end = end
        self.count = count



def read(path):
    countData = []
    binnedData = []
    difficult = []
    mixed = []
    if path[-3:] == ".gz":
        with gzip.open(path, mode='rt') as f:
            for line in f:
                if len(line.strip()) == 0:
                    continue
                elif (line[:2] == "##"):
                    continue
                elif (line[:1] == "#" ):
                    #Generate the matricies for each sample 
                    s = line.strip().split('\t')
                    sampleCount = len(s) - 9
                    #print(s)
                    binnedData = np.zeros(39)
                    for sample in s[9:]:
                        if(sample[0] == '.'):
                            sample = sample[2:]
                        countData.append([sample])
                    continue

                s = line.strip().split('\t')

                #Parse the ID field to detect SNV, Inserts, and Deletes
                variantSize = []
                IDField = s[2].split('-')[0]
                #print(IDField)

                #Check small values
                if (IDField == 'X') or (IDField == 'XX'):
                    variantSize.append(0)
                elif (IDField == 'I'):
                    variantSize.append(1)
                elif (IDField == 'II'):
                    variantSize.append(2)
                elif (IDField == 'D'):
                    variantSize.append(-1)
                elif (IDField == 'DD'):
                    variantSize.append(-2)

                #Check for mixed insertions and deletetions
                elif all(x in IDField for x in ['I', 'D']) or all(x in IDField for x in ['I', 'X']) or all(x in IDField for x in ['X', 'D']) or all(x in IDField for x in ['Y', 'D']) or all(x in IDField for x in ['Y', 'I']):
                    mixed.append(IDField)
                    #print(IDField, tokens)
                    if (IDField == "DX"):
                        variantSize.append(-1)
                        variantSize.append(0)
                    elif (IDField == "DDX"):
                        variantSize.append(-2)
                        variantSize.append(0)
                    elif (IDField == 'IX') or (IDField == 'YX'):
                        variantSize.append(1)
                        variantSize.append(0)
                    elif (IDField == 'IIX') or (IDField == 'YYX'):
                        variantSize.append(2)
                        variantSize.append(0)
                    elif any(x in IDField for x in ['+', '_']):
                        variantSize.append(2000000)
                        difficult.append(IDField)
                    else:
                        tokens = re.finditer("(\d*)([A-Z]+)", IDField)
                        for token in tokens:
                            digit = token.group(1)
                            letters = token.group(2)
                            if (digit == ''):
                                #print(IDField, digit, letters)
                                if (len(letters) == 1):
                                    if (letters == 'I') or (letters == 'Y'):
                                        variantSize.append(1)
                                    elif (letters == 'X') :
                                        variantSize.append(0)
                                    elif (letters == 'D'):
                                        variantSize.append(-1)

                                elif (len(letters) == 2):
                                    if (letters == "II") or (letters[2:] == "YY"):
                                        variantSize.append(2)
                                    elif (letters == "XX"):
                                        variantSize.append(0)
                                    elif (letters == "DD"):
                                        variantSize.append(-2)

                            else:
                                if (letters[0] == 'I') or (letters[0] == 'Y'):
                                    variantSize.append(int(digit))
                                elif (letters[0] == 'X') :
                                    variantSize.append(0)
                                elif (letters[0] == 'D'):
                                    variantSize.append(-int(digit))

                                if (len(letters) == 2):
                                    if (letters[1] == 'I') or (letters[1] == 'Y'):
                                        variantSize.append(1)
                                    elif (letters[1] == 'X') :
                                        variantSize.append(0)
                                    elif (letters[1] == 'D'):
                                        variantSize.append(-1)

                                elif (len(letters) == 3):
                                    if (letters[1:] == "II") or (letters[2:] == "YY"):
                                        variantSize.append(2)
                                    elif (letters[1:] == "XX"):
                                        variantSize.append(0)
                                    elif (letters[1:] == "DD"):
                                        variantSize.append(-2)
                
                #Difficult to parse changes go here.
                elif any(x in IDField for x in ["BND", '+', '_', "DUP"]):
                    variantSize.append(2000000)
                    difficult.append(IDField)

                #Single insertions or deletions of various sizes.
                elif any(x in IDField for x in ['I', 'D', 'Y', 'X']):
                    stringPos = 0
                    for char in IDField:
                        if (char == "I") or (char == 'Y'):
                            variantSize.append(int(IDField[0:stringPos]))
                            break
                        elif (char == "D"):
                            print(IDField)
                            variantSize.append(-int(IDField[0:stringPos]))
                            break
                        elif (char == 'X'):
                            variantSize.append(0)
                            break
                        stringPos += 1      
                
                else:
                    #print("Difficult:", IDField)
                    variantSize.append(2000000)
                    difficult.append(IDField)

                #Update count and binned data
                for sample in np.arange(0, len(s) - 9):
                    GTField = s[9 + sample].split(':')
                    #print(GTField)
                    if (GTField[0] != "0/0"):
                        #print(countData)
                        countData[sample].append(variantSize)
                        for variant in variantSize:
                            print(variant, binFinder.findBin(variant))
                            binnedData[binFinder.findBin(variant)] += 1




                

    return countData, binnedData, difficult, mixed
                

def mergeVcf(path):
    fileNames = []
    for root, dirs, files in os.walk(path):
        for filename in files:
            if(filename[-13:] == ".FINAL.vcf.gz"):
                fileNames.append(root + filename)
    

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
    opts, args = getopt.getopt(argv, "ho:p:m:", ['--path'])
    outputPrefix = ""

    for opt, arg in opts:
        if opt == '-h':
            print("I need help too.")
            sys.exit(1)
        elif opt in ('-p', '--path'):
            path = arg
        elif opt in ('-o', '--output'):
            outputPrefix = arg
    
    countsDict = {}
    binnedDict = {}
    binnedValues = np.zeros(39)
    difficultToParse = []
    mixedParse = []
    
    #Iterate through all ".FINAL.vcf.gz" files and generate basic counts
    print("Reading all vcf files in path:", path)
    for root, dirs, files in os.walk(path):
        for filename in files:
            if(filename[-13:] == ".FINAL.vcf.gz"):
                rawData, binnedData, difficult, mixed = read(root + filename)

                #Merge raw counts with dictionary
                for row in rawData:
                    if row[0] in countsDict:
                        countsDict[row[0]] = np.append(countsDict[row[0]], row[1:])
                    else:
                        countsDict[row[0]] = row[1:]

                #Merge the binned counts with dictionary
                for i in np.arange(0, 39):
                    binnedValues[i] += binnedData[i]

                #Merge counters for difficult or mixed calls.
                if (len(difficult) != 0):
                    difficultToParse.append(difficult)
                if (len(mixed) != 0):
                    mixedParse.append(mixed)

    #bins = np.array([-5000, -2000, -1000, -500, -200, -100, -50, -20, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100, 200, 500, 1000, 2000])

    wCount = csv.writer(open(outputPrefix + "counts.csv", "w"))
    for key, val in countsDict.items():
        wCount.writerow([key,val])

    np.savetxt(outputPrefix + "binned.csv", binnedValues, delimiter=",", fmt="%10.0f")
    np.savetxt(outputPrefix + "mixed.csv", np.asarray(mixedParse), delimiter=",", fmt="%s")
    np.savetxt(outputPrefix + "difficult.csv", np.asarray(difficultToParse), delimiter=",", fmt="%s")

    


if __name__ == '__main__':
      main(sys.argv[1:])
