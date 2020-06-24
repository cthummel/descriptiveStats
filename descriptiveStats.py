import numpy as np
import math, sys, getopt, gzip, os, csv, binFinder, re



def read(path):
    countData = []
    binnedData = []
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
                difficult = []
                mixed = []
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
                    stringPos = 0
                    lastNumberPos = 0
                    mixed.append(IDField)
                    tokens = re.findall("(\d*)([A-Z]+)", IDField)
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
                    else:
                        tokens = re.finditer("(\d*)([A-Z]+)", IDField)
                        for token in tokens:
                            digit = token.group(1)
                            letters = token.group(2)
                            if (digit == ''):
                                print(IDField, digit, letters)

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
                

                #Single insertions or deletions of various sizes.
                elif any(x in IDField for x in ['I', 'D', 'Y', 'X']):
                    stringPos = 0
                    for char in IDField:
                        if (char == "I") or (char == 'Y'):
                            variantSize.append(int(IDField[0:stringPos]))
                            break
                        elif (char == "D") and not IDField.find("BND"):
                            variantSize.append(-int(IDField[0:stringPos]))
                            break
                        elif (char == 'X'):
                            variantSize.append(0)
                            break
                        stringPos += 1
                #Difficult to parse changes go here.
                elif any(x in IDField for x in ["BND"]):
                    variantSize.append(2000000)
                    difficult.append(IDField)
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
                            binnedData[binFinder.findBin(variant)] += 1




                

    return countData, binnedData, difficult, mixed
                


def main(argv):
    opts, args = getopt.getopt(argv, "ho:p:", ['--path'])
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
                #print(filename)
                #with open(os.path.join(root, filename), 'r') as f:
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
                # for row in binnedData:
                #     if row[0] in binnedDict:
                #         for value in row[1:]:
                #             binnedDict[row[0]] += value
                #     else:
                #         binnedDict[row[0]] = row[1:]
                if (len(difficult) != 0):
                    difficultToParse.append(difficult)
                if (len(mixed) != 0):
                    mixedParse.append(mixed)

    #bins = np.array([-5000, -2000, -1000, -500, -200, -100, -50, -20, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100, 200, 500, 1000, 2000])

    wCount = csv.writer(open(outputPrefix + "counts.csv", "w"))
    for key, val in countsDict.items():
        wCount.writerow([key, val])

    np.savetxt(outputPrefix + "binned.csv", binnedValues, delimiter=",", fmt="%10.0f")
    np.savetxt(outputPrefix + "mixed.csv", np.asarray(mixedParse), delimiter=",", fmt="%s")
    np.savetxt(outputPrefix + "difficult.csv", np.asarray(difficultToParse), delimiter=",", fmt="%s")

    


if __name__ == '__main__':
      main(sys.argv[1:])
