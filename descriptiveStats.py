import numpy as np
import math, sys, getopt, gzip, os, csv




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
                    print(s)
                    binnedData = np.zeros((len(s) - 9, 39))
                    for sample in s[9:]:
                        if(sample[0] == '.'):
                            sample = sample[2:]
                        countData.append([sample])
                    continue

                s = line.strip().split('\t')

                #Parse the ID field to detect SNV, Inserts, and Deletes
                variantSize = 0
                IDField = s[2].split('-')[0]
                print(IDField)

                #Check small values
                if (IDField == 'X'):
                    variantSize = 0
                elif (IDField == 'I'):
                    variantSize = 1
                elif (IDField == 'II'):
                    variantSize = 2
                elif (IDField == 'D'):
                    variantSize = -1
                elif (IDField == 'DD'):
                    variantSize = -2
                #Check for mixed insertions and deletetions
                elif all(x in IDField for x in ['I', 'D']):
                    print("Mixed:", IDField)
                    variantSize = 1000000
                #     if 'I' in IDField:

                #     else:

                #Single insertions or deletions of various sizes.
                elif any(x in IDField for x in ['I', 'D']):
                    stringPos = 0
                    for char in IDField:
                        if (char == "I"):
                            variantSize = int(IDField[0:stringPos])
                            break
                        elif (char == "D"):
                            variantSize = -int(IDField[0:stringPos])
                            break
                        stringPos += 1
                #Difficult to parse changes go here.
                else:
                    print("Difficult:", IDField)
                    variantSize = 2000000


                

                

                #Update count and binned data
                for sample in np.arange(0, len(s) - 9):
                    GTField = s[9 + sample].split(':')
                    if (GTField[0] != "0/0"):
                        print(countData)
                        countData[sample].append(variantSize)




                

    return countData, binnedData
                


def main(argv):
    opts, args = getopt.getopt(argv, "hop:", ['--path'])
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
    
    #Iterate through all ".FINAL.vcf.gz" files and generate basic counts
    for root, dirs, files in os.walk(path):
        for filename in files:
            if(filename[-13:] == ".FINAL.vcf.gz"):
                print(filename)
                #with open(os.path.join(root, filename), 'r') as f:
                rawData, binnedData = read(root + filename)

                #Merge raw counts with dictionary
                for row in rawData:
                    if row[0] in countsDict:
                        countsDict[row[0]] = np.append(countsDict[row[0]], row[1:])
                    else:
                        countsDict[row[0]] = row[1:]

                #Merge the binned counts with dictionary
                for row in binnedData:
                    if row[0] in binnedDict:
                        for value in row[1:]:
                            binnedDict[row[0]] += value
                    else:
                        binnedDict[row[0]] = row[1:]


    wCount = csv.writer(open(outputPrefix + "counts.csv", "w"))
    for key, val in countsDict.items():
        wCount.writerow([key, val])

    wBinned = csv.writer(open(outputPrefix + "binned.csv", "w"))
    for key, val in binnedDict.items():
        wBinned.writerow([key, val])

    


if __name__ == '__main__':
      main(sys.argv[1:])
