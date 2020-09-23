import sys, getopt, gzip, csv, fileinput, os, statistics
import numpy as np

class megabaseInfo:
    def __init__(self, start, end, count, insertion, deletion, snv, fatherAge, motherAge, varPos):
        self.start = start
        self.end = end
        self.count = count
        self.insertion = insertion
        self.deletion = deletion
        self.snv = snv
        self.fatherAge = fatherAge
        self.motherAge = motherAge
        self.variantPosition = varPos

class geneInfo:
    def __init__(self, name, transcriptStart, transcriptEnd, count, insertion, deletion, snv, fatherAge, motherAge, varPos):
        self.name = name
        self.transcriptStart = transcriptStart
        self.transcriptEnd = transcriptEnd
        # self.exonStart = exonStart
        # self.exonEnd = exonEnd
        self.count = count
        self.insertion = insertion
        self.deletion = deletion
        self.snv = snv
        self.fatherAge = fatherAge
        self.motherAge = motherAge
        self.variantPosition = varPos


class familyInfo:
    def __init__(self, familyID, probandID, siblingID, probandGender, siblingGender, siblingMotherAge, siblingFatherAge, probandMotherAge, probandFatherAge):
        self.familyID = familyID
        self.probandID = probandID
        self.siblingID = siblingID
        self.probandGender = probandGender
        self.siblingGender = siblingGender
        self.siblingFatherAge = siblingFatherAge
        self.siblingMotherAge = siblingMotherAge
        self.probandFatherAge = probandFatherAge
        self.probandMotherAge = probandMotherAge
        
def ageMean(data):
    result = 0
    count = 0
    for x in data:
        if x == "NA":
            continue
        else:
            result += x
            count += 1

    if count == 0:
        return "NA"
    else:
        return result * 1.0 / count

def pairSiblings(famFile, sampleFile):
    families = []

    with open(sampleFile, mode='rt') as f:
        header = f.readline()
        for line in f:
            s = line.strip().split("\t")
            families.append(familyInfo(s[0], s[3], s[4], "?", "?", s[16], s[17], s[18], s[19]))

    familyIndex = 0
    currentFamilyId = 0

    with open(famFile, mode='rt') as f:
        header = f.readline()
        for line in f:
            if familyIndex == len(families):
                break
            s = line.strip().split()
            ind = s[1].split(".")[1]
            if int(s[0]) != int(families[familyIndex].familyID):
                while int(s[0]) > int(families[familyIndex].familyID):
                    familyIndex += 1
                    if familyIndex == len(families):
                        break
                if familyIndex == len(families):
                        break
                if (int(s[0]) < int(families[familyIndex].familyID)):
                    continue
            else:
                if ind == "p1":
                    families[familyIndex].probandGender = s[2]
                elif ind == "s1":
                    families[familyIndex].siblingGender = s[2]
                
                # if families[familyIndex].siblingGender != "?" and families[familyIndex].probandGender != "?":
                #     print(families[familyIndex].siblingGender, families[familyIndex].probandGender)

    return families


def appendChromInfoDict(data, megaBaseStart, megaBaseSize, chrom, variantCount, insert, Del, snv, fatherAge, motherAge, varPos):
    if (fatherAge == [] or motherAge == []):
        data[chrom].append(megabaseInfo(megaBaseStart, megaBaseStart + megaBaseSize, variantCount, insert, Del, snv, [], [], []))
    else:
        data[chrom].append(megabaseInfo(megaBaseStart, megaBaseStart + megaBaseSize, variantCount, insert, Del, snv, [int(fatherAge)], [int(motherAge)], [varPos]))

def updateChromInfoDict(data, index, chrom, variantCount, insert, Del, snv, fatherAge, motherAge, varPos, updateNext):
    data[chrom][index].count += variantCount
    data[chrom][index].insertion += insert
    data[chrom][index].deletion += Del
    data[chrom][index].snv += snv
    data[chrom][index].variantPosition.append(varPos)
    if fatherAge != "":
        data[chrom][index].fatherAge.append(int(fatherAge))
    else:
        data[chrom][index].fatherAge.append("NA")
    if motherAge != "":
        data[chrom][index].motherAge.append(int(motherAge))
    else:
        data[chrom][index].motherAge.append("NA")
    if updateNext:
        data[chrom][index + 1].count += variantCount
        data[chrom][index + 1].insertion += insert
        data[chrom][index + 1].deletion += Del
        data[chrom][index + 1].snv += snv
        data[chrom][index + 1].variantPosition.append(varPos)
        if fatherAge != "":
            data[chrom][index + 1].fatherAge.append(int(fatherAge))
        else:
            data[chrom][index + 1].fatherAge.append("NA")
        if motherAge != "":
            data[chrom][index + 1].motherAge.append(int(motherAge))
        else:
            data[chrom][index + 1].motherAge.append("NA")


def geneCountMergeFamily(file, outputPrefix, familyData):
    probandDataSet = True
    # chromInfoDict = {}
    # chromInfoDictMM = {}
    # chromInfoDictMF = {}
    # chromInfoDictFM = {}
    # chromInfoDictFF = {}
    #[MM, MF, FM, FF, male, female]
    result = [{}, {}, {}, {}, {}, {}]

    #generate bins
    with gzip.open("gencode.v34.annotation.gff3.gz", mode='rt') as f:
        for line in f:
            
            s = line.strip().split("\t")
            if s[0][0] == "#":
                continue

            infoField = s[8].strip().split(";")
             
            if s[2] == "gene" and infoField[2][10:] == "protein_coding": #gene_type=
                geneName = infoField[3][10:]
                for i in np.arange(0, len(result)):
                    if (s[0] not in result[i].keys()):
                        result[i][s[0]] = [geneInfo(geneName, int(s[3]), int(s[4]), 0, 0, 0, 0, [], [], [])]
                    else:
                        result[i][s[0]].append(geneInfo(geneName, int(s[3]), int(s[4]), 0, 0, 0, 0, [], [], []))


    for root, dirs, files in os.walk(file):
        temproot = root.strip().split("/")
        #print("root", temproot)
        if temproot[len(temproot) - 2] == "s1":
            probandDataSet = False
        for filename in files:
            if(filename[-13:] == ".FINAL.vcf.gz"):
                with gzip.open(root + filename, mode='rt') as f:
                    print("Scanning variants from file:", filename)
                    currentChrom = ""
                    currentMegaBaseIndex = 0
                    currentDataSet = 0
                    currentFatherAge = 0
                    currentMotherAge = 0
                    gender = 4
                    for line in f:
                        if len(line.strip()) == 0:
                            continue
                        elif (line[:2] == "##"):
                            continue
                        elif (line[:1] == "#"):
                            #parse which sample the current file has.
                            sample = line.strip().split('\t')[9][2:]
                            #print("Sample:", sample)
                            if probandDataSet:
                                for x in familyData:
                                    if x.probandID == sample:
                                        currentFatherAge = x.probandFatherAge
                                        currentMotherAge = x.probandMotherAge
                                        if (x.probandGender == "male") and (x.siblingGender == "male"):
                                            currentDataSet = 0
                                        elif (x.probandGender == "male") and (x.siblingGender == "female"):
                                            currentDataSet = 1
                                        elif (x.probandGender == "female") and (x.siblingGender == "male"):
                                            currentDataSet = 2
                                            gender = 5
                                        elif (x.probandGender == "female") and (x.siblingGender == "female"):
                                            currentDataSet = 3
                                            gender = 5
                                        #print(currentDataSet, x.probandGender, x.siblingGender)
                                        break            
                            else:
                                for x in familyData:
                                    if x.siblingID == sample:
                                        currentFatherAge = x.siblingFatherAge
                                        currentMotherAge = x.siblingMotherAge
                                        if (x.probandGender == "male") and (x.siblingGender == "male"):
                                            currentDataSet = 0
                                        elif (x.probandGender == "male") and (x.siblingGender == "female"):
                                            currentDataSet = 1
                                        elif (x.probandGender == "female") and (x.siblingGender == "male"):
                                            currentDataSet = 2
                                            gender = 5
                                        elif (x.probandGender == "female") and (x.siblingGender == "female"):
                                            currentDataSet = 3
                                            gender = 5
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
                            currentMegaBaseEnd = result[currentDataSet][s[0]][0].transcriptEnd
                            currentChrom = s[0]

                        for gene in result[currentDataSet][currentChrom]:
                            if gene.transcriptStart <= int(s[1]) and gene.transcriptEnd >= int(s[1]):
                                gene.count += variantCount
                                gene.insertion += insert
                                gene.deletion += Del
                                gene.snv += snv
                                gene.variantPosition.append(int(s[1]))
                                if currentFatherAge != "":
                                    gene.fatherAge.append(int(currentFatherAge))
                                else:
                                    gene.fatherAge.append("NA")
                                if currentMotherAge != "":
                                    gene.motherAge.append(int(currentMotherAge))
                                else:
                                    gene.motherAge.append("NA")
                            elif gene.transcriptStart > int(s[1]):
                                break
                        
                        #for gender combined dataset
                        for gene in result[gender][currentChrom]:
                            if gene.transcriptStart <= int(s[1]) and gene.transcriptEnd >= int(s[1]):
                                gene.count += variantCount
                                gene.insertion += insert
                                gene.deletion += Del
                                gene.snv += snv
                                gene.variantPosition.append(int(s[1]))
                                if currentFatherAge != "":
                                    gene.fatherAge.append(int(currentFatherAge))
                                else:
                                    gene.fatherAge.append("NA")
                                if currentMotherAge != "":
                                    gene.motherAge.append(int(currentMotherAge))
                                else:
                                    gene.motherAge.append("NA")
                            elif gene.transcriptStart > int(s[1]):
                                break
                            
                        #updateChromInfoDict(result[currentDataSet], currentMegaBaseIndex, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, False)
                        

    outputIndex = 0
    typePrefix = ["MM", "MF", "FM", "FF", "male", "female"]

    # maleCounts = csv.writer(open(outputPrefix + "male" + ".geneCounts.csv", "w"))
    # maleAgeCounts = csv.writer(open(outputPrefix + "male" + ".geneAgeVector.csv", "w"), delimiter="\t")
    # femaleCounts = csv.writer(open(outputPrefix + "female" + ".geneCounts.csv", "w"))
    # femaleAgeCounts = csv.writer(open(outputPrefix + "female" + ".geneAgeVector.csv", "w"), delimiter="\t")

    # maleCounts.writerow(["Chrom", "Start", "End", "Count", "Insertions", "Deletions", "SNV", "MeanFatherAge", "MeanMotherAge", "Gene"])
    # maleAgeCounts.writerow(["Chrom", "Gene", "Start", "End", "FatherAge", "MotherAge", "VariantPosition"])
    # femaleCounts.writerow(["Chrom", "Start", "End", "Count", "Insertions", "Deletions", "SNV", "MeanFatherAge", "MeanMotherAge", "Gene"])
    # femaleAgeCounts.writerow(["Chrom", "Gene", "Start", "End", "FatherAge", "MotherAge", "VariantPosition"])

    for x in result:
        wCounts = csv.writer(open(outputPrefix + typePrefix[outputIndex] + ".geneCounts.csv", "w"))
        fAgeCounts = csv.writer(open(outputPrefix + typePrefix[outputIndex] + ".geneAgeVector.csv", "w"), delimiter="\t")
        fAgeCounts.writerow(["Chrom", "Gene", "Start", "End", "FatherAge", "MotherAge", "VariantPosition"])
        wCounts.writerow(["Chrom", "Start", "End", "Count", "Insertions", "Deletions", "SNV", "MeanFatherAge", "MeanMotherAge", "Gene"])
        for key in sorted(x.keys()):
            for val in x[key]:
                fAgeCounts.writerow([key, val.name, val.transcriptStart, val.transcriptEnd, val.fatherAge, val.motherAge, val.variantPosition])
                # if outputIndex < 2:
                #     maleAgeCounts.writerow([key, val.name, val.transcriptStart, val.transcriptEnd, val.fatherAge, val.motherAge, val.variantPosition])
                # else:
                #     femaleAgeCounts.writerow([key, val.name, val.transcriptStart, val.transcriptEnd, val.fatherAge, val.motherAge, val.variantPosition])

                if val.fatherAge == [] or val.motherAge == []:
                    wCounts.writerow([key, val.transcriptStart, val.transcriptEnd, val.count, val.insertion, val.deletion, val.snv, "NA", "NA", val.name])
                    # if outputIndex < 2:
                    #     maleCounts.writerow([key, val.transcriptStart, val.transcriptEnd, val.count, val.insertion, val.deletion, val.snv, "NA", "NA", val.name])
                    # else:
                    #     femaleCounts.writerow([key, val.transcriptStart, val.transcriptEnd, val.count, val.insertion, val.deletion, val.snv, "NA", "NA", val.name])
                else:
                    wCounts.writerow([key, val.transcriptStart, val.transcriptEnd, val.count, val.insertion, val.deletion, val.snv, ageMean(val.fatherAge), ageMean(val.motherAge), val.name])
                    # if outputIndex < 2:
                    #     maleCounts.writerow([key, val.transcriptStart, val.transcriptEnd, val.count, val.insertion, val.deletion, val.snv, ageMean(val.fatherAge), ageMean(val.motherAge), val.name])
                    # else:
                    #     femaleCounts.writerow([key, val.transcriptStart, val.transcriptEnd, val.count, val.insertion, val.deletion, val.snv, ageMean(val.fatherAge), ageMean(val.motherAge), val.name])

        outputIndex += 1

    return result

    
def megabaseCountMergeFamily(file, overlap, binsize, outputPrefix, familyData):
    megabaseSize = binsize
    probandDataSet = True
    chromInfoDict = {}
    chromInfoDictMM = {}
    chromInfoDictMF = {}
    chromInfoDictFM = {}
    chromInfoDictFF = {}
    #[MM, MF, FM, FF, male, female]
    result = [{}, {}, {}, {}, {}, {}]

    for root, dirs, files in os.walk(file):
        temproot = root.strip().split("/")
        #print("root", temproot)
        if temproot[len(temproot) - 2] == "s1":
            probandDataSet = False
        for filename in files:
            if(filename[-13:] == ".FINAL.vcf.gz"):
                with gzip.open(root + filename, mode='rt') as f:
                    print("Scanning variants from file:", filename)
                    currentChrom = ""
                    currentMegaBaseIndex = 0
                    genderCurrentMegaBaseIndex = 0
                    currentDataSet = 0
                    currentFatherAge = 0
                    currentMotherAge = 0
                    gender = 4
                    for line in f:
                        if len(line.strip()) == 0:
                            continue
                        elif (line[:2] == "##"):
                            continue
                        elif (line[:1] == "#"):
                            #parse which sample the current file has.
                            sample = line.strip().split('\t')[9][2:]
                            #print("Sample:", sample)
                            if probandDataSet:
                                for x in familyData:
                                    if x.probandID == sample:
                                        currentFatherAge = x.probandFatherAge
                                        currentMotherAge = x.probandMotherAge
                                        if (x.probandGender == "male") and (x.siblingGender == "male"):
                                            currentDataSet = 0
                                        elif (x.probandGender == "male") and (x.siblingGender == "female"):
                                            currentDataSet = 1
                                        elif (x.probandGender == "female") and (x.siblingGender == "male"):
                                            currentDataSet = 2
                                            gender = 5
                                        elif (x.probandGender == "female") and (x.siblingGender == "female"):
                                            currentDataSet = 3
                                            gender = 5
                                        #print(currentDataSet, x.probandGender, x.siblingGender)
                                        break            
                            else:
                                for x in familyData:
                                    if x.siblingID == sample:
                                        currentFatherAge = x.siblingFatherAge
                                        currentMotherAge = x.siblingMotherAge
                                        if (x.probandGender == "male") and (x.siblingGender == "male"):
                                            currentDataSet = 0
                                        elif (x.probandGender == "male") and (x.siblingGender == "female"):
                                            currentDataSet = 1
                                        elif (x.probandGender == "female") and (x.siblingGender == "male"):
                                            currentDataSet = 2
                                            gender = 5
                                        elif (x.probandGender == "female") and (x.siblingGender == "female"):
                                            currentDataSet = 3
                                            gender = 5
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
                                result[currentDataSet][s[0]] = [megabaseInfo(0, 0 + megabaseSize, 0, 0, 0, 0, [], [], [])]    
                                currentMegaBaseEnd = 0 + megabaseSize
                            
                            if (s[0] not in result[gender]):
                                result[gender][s[0]] = [megabaseInfo(0, 0 + megabaseSize, 0, 0, 0, 0, [], [], [])]    
                                genderCurrentMegaBaseEnd = 0 + megabaseSize
                            
                            currentMegaBaseIndex = 0
                            currentMegaBaseEnd = result[gender][s[0]][0].end
                            currentChrom = s[0]

                            genderCurrentMegaBaseIndex = 0
                            genderCurrentMegaBaseEnd = result[gender][s[0]][0].end
                            genderCurrentChrom = s[0]

                        while (int(s[1]) > currentMegaBaseEnd):
                            #If the megabase isnt initialized yet
                            if (currentMegaBaseIndex + 1 == len(result[currentDataSet][s[0]])):
                                newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
                                #chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, 0, 0, 0, 0))
                                appendChromInfoDict(result[currentDataSet], newMegaBaseStart, megabaseSize, s[0], 0, 0, 0, 0, [], [], [])
                                currentMegaBaseEnd = newMegaBaseStart + megabaseSize
                                currentMegaBaseIndex += 1
                            else:
                                currentMegaBaseIndex += 1
                                currentMegaBaseEnd = result[currentDataSet][s[0]][currentMegaBaseIndex].end

                        #If we find a variant that falls in the overlap for two megabases.
                        if (int(s[1]) <= currentMegaBaseEnd and int(s[1]) > currentMegaBaseEnd - (megabaseSize * overlap)):
                            if (currentMegaBaseIndex == len(result[currentDataSet][s[0]]) - 1):
                                newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
                                appendChromInfoDict(result[currentDataSet], newMegaBaseStart, megabaseSize, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, int(s[1]))
                                #chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, variantCount, insert, Del, snv))
                                updateChromInfoDict(result[currentDataSet], currentMegaBaseIndex, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, int(s[1]), False)
                            else:
                                updateChromInfoDict(result[currentDataSet], currentMegaBaseIndex, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, int(s[1]), True)
                        else:
                            updateChromInfoDict(result[currentDataSet], currentMegaBaseIndex, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, int(s[1]), False)
                            

                        #Do it again to merge properly with whole gender data set   
                        while (int(s[1]) > genderCurrentMegaBaseEnd):
                            #If the megabase isnt initialized yet
                            if (genderCurrentMegaBaseIndex + 1 == len(result[gender][s[0]])):
                                newMegaBaseStart = genderCurrentMegaBaseEnd - (overlap * megabaseSize)
                                #chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, 0, 0, 0, 0))
                                appendChromInfoDict(result[gender], newMegaBaseStart, megabaseSize, s[0], 0, 0, 0, 0, [], [], [])
                                genderCurrentMegaBaseEnd = newMegaBaseStart + megabaseSize
                                genderCurrentMegaBaseIndex += 1
                            else:
                                genderCurrentMegaBaseIndex += 1
                                genderCurrentMegaBaseEnd = result[gender][s[0]][currentMegaBaseIndex].end
                                

                        #If we find a variant that falls in the overlap for two megabases.
                        if (int(s[1]) <= genderCurrentMegaBaseEnd and int(s[1]) > genderCurrentMegaBaseEnd - (megabaseSize * overlap)):
                            if (currentMegaBaseIndex == len(result[gender][s[0]]) - 1):
                                newMegaBaseStart = genderCurrentMegaBaseEnd - (overlap * megabaseSize)
                                appendChromInfoDict(result[gender], newMegaBaseStart, megabaseSize, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, int(s[1]))
                                #chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, variantCount, insert, Del, snv))
                                updateChromInfoDict(result[gender], genderCurrentMegaBaseIndex, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, int(s[1]), False)
                            else:
                                updateChromInfoDict(result[gender], genderCurrentMegaBaseIndex, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, int(s[1]), True)
                        else:
                            updateChromInfoDict(result[gender], genderCurrentMegaBaseIndex, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, int(s[1]), False)
                            
                        
                        #print(chrom, int(s[1]), chromInfoDict[chrom][currentMegaBaseIndex].end - megabaseSize, chromInfoDict[chrom][currentMegaBaseIndex].end, chromInfoDict[chrom][currentMegaBaseIndex].count)

    # maleCounts = csv.writer(open(outputPrefix + "male" + ".megaBaseCounts.csv", "w"))
    # maleAgeCounts = csv.writer(open(outputPrefix + "male" + ".megaBaseAgeVector.csv", "w"), delimiter="\t")
    # femaleCounts = csv.writer(open(outputPrefix + "female" + ".megaBaseCounts.csv", "w"))
    # femaleAgeCounts = csv.writer(open(outputPrefix + "female" + ".megaBaseAgeVector.csv", "w"), delimiter="\t")

    # maleCounts.writerow(["Chrom", "Start", "End", "Count", "Insertions", "Deletions", "SNV", "MeanFatherAge", "MeanMotherAge", "Gene"])
    # maleAgeCounts.writerow(["Chrom", "Gene", "Start", "End", "FatherAge", "MotherAge", "VariantPosition"])
    # femaleCounts.writerow(["Chrom", "Start", "End", "Count", "Insertions", "Deletions", "SNV", "MeanFatherAge", "MeanMotherAge", "Gene"])
    # femaleAgeCounts.writerow(["Chrom", "Gene", "Start", "End", "FatherAge", "MotherAge", "VariantPosition"])

    outputIndex = 0
    typePrefix = ["MM", "MF", "FM", "FF", "male", "female"]
    for x in result:
        wCounts = csv.writer(open(outputPrefix + typePrefix[outputIndex] + ".megaBaseCounts.csv", "w"))
        mAgeCounts = csv.writer(open(outputPrefix + typePrefix[outputIndex] + ".megaBaseAgeVector.csv", "w"), delimiter="\t")
        mAgeCounts.writerow(["Chrom", "Gene", "Start", "End", "FatherAge", "MotherAge", "VariantPosition"])
        wCounts.writerow(["Chrom", "Start", "End", "Count", "Insertions", "Deletions", "SNV", "MeanFatherAge", "MeanMotherAge", "Gene"])
        for key in sorted(x.keys()):
            for val in x[key]:
                #print(val.fatherAge)
                mAgeCounts.writerow([key, "NA", val.start, val.end, val.fatherAge, val.motherAge, val.variantPosition])
                if val.fatherAge == [] or val.motherAge == []:
                    wCounts.writerow([key, int(val.start), int(val.end), val.count, val.insertion, val.deletion, val.snv, "NA", "NA", "NA"])
                else:
                    wCounts.writerow([key, int(val.start), int(val.end), val.count, val.insertion, val.deletion, val.snv, ageMean(val.fatherAge), ageMean(val.motherAge), "NA"])
                    
        outputIndex += 1

    return result




# def megabaseCountMerge(file, overlap, binsize, outputPrefix):
#     megabaseSize = binsize
#     chromInfoDict = {}
#     for root, dirs, files in os.walk(file):
#         for filename in files:
#             if(filename[-13:] == ".FINAL.vcf.gz"):
#                 with gzip.open(root + filename, mode='rt') as f:
#                     print("Scanning variants from file:", filename)
#                     currentChrom = ""
#                     currentMegaBaseIndex = 0
#                     for line in f:
#                         if len(line.strip()) == 0:
#                             continue
#                         elif (line[:2] == "##"):
#                             continue
#                         elif (line[:1] == "#"):
#                             continue
#                         s = line.strip().split('\t')
#                         IDField = s[2].split('-')[0]

#                         insert = 0
#                         Del = 0
#                         snv = 0

#                         #Parse if the variant is an insertion, deletion or both
#                         if any(x in IDField for x in ["BND", '+', '_', "DUP"]):
#                             continue
#                         else:
#                             if ("I" in IDField):
#                                 insert += 1
#                             if ("Y" in IDField):
#                                 insert += 1
#                             if ("D" in IDField):
#                                 Del += 1
#                             if ("X" in IDField):
#                                 snv += 1

#                         variantCount = insert + Del + snv

#                         #Place the variant in the right megabase bin
#                         #New chromosome means we add a new key to the dictionary and append a new megabase counter.
#                         ##s[0] is the chromosome, s[1] is the position
#                         if (len(s[0]) > 5):
#                             continue
#                         if (s[0] != currentChrom):
#                             if (s[0] not in chromInfoDict):
#                                 chromInfoDict[s[0]] = [megabaseInfo(0, 0 + megabaseSize, 0, 0, 0, 0, 0, 0)]    
#                                 currentMegaBaseEnd = 0 + megabaseSize
                            
#                             currentMegaBaseIndex = 0
#                             currentMegaBaseEnd = chromInfoDict[s[0]][0].end
#                             currentChrom = s[0]

#                         while (int(s[1]) > currentMegaBaseEnd):
#                             #If the megabase isnt initialized yet
#                             if (currentMegaBaseIndex + 1 == len(chromInfoDict[s[0]])):
#                                 newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
#                                 chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, 0, 0, 0, 0, 0, 0))
#                                 currentMegaBaseEnd = newMegaBaseStart + megabaseSize
#                                 currentMegaBaseIndex += 1
#                             else:
#                                 currentMegaBaseIndex += 1
#                                 currentMegaBaseEnd = chromInfoDict[s[0]][currentMegaBaseIndex].end
                                

#                         #If we find a variant that falls in the overlap for two megabases.
#                         if (int(s[1]) <= currentMegaBaseEnd and int(s[1]) > currentMegaBaseEnd - (megabaseSize * overlap)):
#                             if (currentMegaBaseIndex == len(chromInfoDict[s[0]]) - 1):
#                                 newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
#                                 chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, variantCount, insert, Del, snv, 0, 0))
#                                 chromInfoDict[s[0]][currentMegaBaseIndex].count += variantCount
#                                 chromInfoDict[s[0]][currentMegaBaseIndex].insertion += insert
#                                 chromInfoDict[s[0]][currentMegaBaseIndex].deletion += Del
#                                 chromInfoDict[s[0]][currentMegaBaseIndex].snv += snv
#                             else:
#                                 chromInfoDict[s[0]][currentMegaBaseIndex].count += variantCount
#                                 chromInfoDict[s[0]][currentMegaBaseIndex].insertion += insert
#                                 chromInfoDict[s[0]][currentMegaBaseIndex].deletion += Del
#                                 chromInfoDict[s[0]][currentMegaBaseIndex].snv += snv
#                                 chromInfoDict[s[0]][currentMegaBaseIndex + 1].count += variantCount
#                                 chromInfoDict[s[0]][currentMegaBaseIndex + 1].insertion += insert
#                                 chromInfoDict[s[0]][currentMegaBaseIndex + 1].deletion += Del
#                                 chromInfoDict[s[0]][currentMegaBaseIndex + 1].snv += snv
#                         else:
#                             chromInfoDict[s[0]][currentMegaBaseIndex].count += variantCount
#                             chromInfoDict[s[0]][currentMegaBaseIndex].insertion += insert
#                             chromInfoDict[s[0]][currentMegaBaseIndex].deletion += Del
#                             chromInfoDict[s[0]][currentMegaBaseIndex].snv += snv
                            
#                         #print(chrom, int(s[1]), chromInfoDict[chrom][currentMegaBaseIndex].end - megabaseSize, chromInfoDict[chrom][currentMegaBaseIndex].end, chromInfoDict[chrom][currentMegaBaseIndex].count)

#     wCounts = csv.writer(open(outputPrefix + "megaBaseCounts.csv", "w"))
#     wCounts.writerow(["Chrom", "Start", "End", "Count", "Insertions", "Deletions", "SNV"])
#     for key in chromInfoDict.keys():
#         for val in chromInfoDict[key]:
#             wCounts.writerow([key, val.start, val.end, val.count, val.insertion, val.deletion, val.snv])

#     return chromInfoDict

# def megabaseCount(file, overlap, binsize, outputPrefix):
#     megabaseSize = 1000000
#     chromInfoDict = {}
#     currentChrom = ""
#     currentMegaBaseIndex = 0
#     for line in fileinput.input(file):
#         if len(line.strip()) == 0:
#             continue
#         elif (line[:2] == "##"):
#             continue
#         elif (line[:1] == "#"):
#             continue
#         s = line.strip().split('\t')
#         #New chromosome means we add a new key to the dictionary and append a new megabase counter.
#         ##s[0] is the chromosome, s[1] is the position
#         chrom = s[0].split('_')[0]
#         IDField = s[2].split('-')[0]

#         insert = 0
#         Del = 0
#         snv = 0

#         #Parse if the variant is an insertion, deletion or both
#         if any(x in IDField for x in ["BND", '+', '_', "DUP"]):
#             continue
#         else:
#             if ("I" in IDField):
#                 insert += 1
#             if ("Y" in IDField):
#                 insert += 1
#             if ("D" in IDField):
#                 Del += 1
#             if ("X" in IDField):
#                 snv += 1

#         variantCount = insert + Del + snv

#         #Place the variant in the right megabase bin
#         if (chrom != currentChrom):
#             print("Scanning variants in:", chrom)
#             if (chrom not in chromInfoDict):
#                 chromInfoDict[s[0]] = [megabaseInfo(0, 0 + megabaseSize, 0, 0, 0, 0, 0, 0)]    
#                 currentMegaBaseEnd = 0 + megabaseSize
            
#             currentMegaBaseEnd = chromInfoDict[s[0]][0].end
#             currentMegaBaseIndex = 0
#             currentChrom = chrom
            

#         while (int(s[1]) > currentMegaBaseEnd):
#             #If the megabase isnt initialized yet
#             if (currentMegaBaseIndex + 1 == len(chromInfoDict[s[0]]) ):
#                 newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
#                 chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, 0, 0, 0 ,0, 0, 0))
#                 currentMegaBaseEnd = newMegaBaseStart + megabaseSize
#                 currentMegaBaseIndex += 1
#             else:
#                 currentMegaBaseIndex += 1
#                 currentMegaBaseEnd = chromInfoDict[s[0]][currentMegaBaseIndex].end
                

#         #If we find a variant that falls in the overlap for two megabases.
#         if (int(s[1]) <= currentMegaBaseEnd and int(s[1]) > currentMegaBaseEnd - (megabaseSize * overlap)):
#             if (currentMegaBaseIndex == len(chromInfoDict[s[0]]) - 1):
#                 newMegaBaseStart = currentMegaBaseEnd - (overlap * megabaseSize)
#                 chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, variantCount, insert, Del, snv, 0, 0))
#                 chromInfoDict[s[0]][currentMegaBaseIndex].count += variantCount
#                 chromInfoDict[s[0]][currentMegaBaseIndex].insertion += insert
#                 chromInfoDict[s[0]][currentMegaBaseIndex].deletion += Del
#                 chromInfoDict[s[0]][currentMegaBaseIndex].snv += snv
#             else:
#                 chromInfoDict[s[0]][currentMegaBaseIndex].count += variantCount
#                 chromInfoDict[s[0]][currentMegaBaseIndex].insertion += insert
#                 chromInfoDict[s[0]][currentMegaBaseIndex].deletion += Del
#                 chromInfoDict[s[0]][currentMegaBaseIndex].snv += snv
#                 chromInfoDict[s[0]][currentMegaBaseIndex + 1].count += variantCount
#                 chromInfoDict[s[0]][currentMegaBaseIndex + 1].insertion += insert
#                 chromInfoDict[s[0]][currentMegaBaseIndex + 1].deletion += Del
#                 chromInfoDict[s[0]][currentMegaBaseIndex + 1].snv += snv
#         else:
#             chromInfoDict[s[0]][currentMegaBaseIndex].count += variantCount
#             chromInfoDict[s[0]][currentMegaBaseIndex].insertion += insert
#             chromInfoDict[s[0]][currentMegaBaseIndex].deletion += Del
#             chromInfoDict[s[0]][currentMegaBaseIndex].snv += snv

#     wCounts = csv.writer(open(outputPrefix + "megaBaseCounts.csv", "w"))
#     wCounts.writerow(["Chrom", "Start", "End", "Count", "Insertions", "Deletions", "SNV"])
#     for key in chromInfoDict.keys():
#         for val in chromInfoDict[key]:
#             wCounts.writerow([key, val.start, val.end, val.count, val.insertion, val.deletion])

#     return chromInfoDict


def main(argv):
    opts, args = getopt.getopt(argv, "h p:", ['merge=', 'output=', 'overlap=', 'binsize=', 'famFile=', 'sampleFile='])
    outputPrefix = ""
    path = '-'
    merge = False
    geneMode = False
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
            if arg == "gene":
                geneMode = True
            else:
                binsize = int(arg)
        elif opt == '--famFile':
            famFile = arg
        elif opt == '--sampleFile':
            sampleFile = arg
        elif opt in ('-h'):
            print("Use", ['--merge', '--output', '--overlap'], "to adjust parameters")
            sys.exit(0)

    familyData = pairSiblings(famFile, sampleFile)

    famAge = csv.writer(open(outputPrefix + "famAge.csv", "w"))
    famAge.writerow(["famID", "probandGender", "siblingGender", "fatherAge"])
    for val in familyData:
        famAge.writerow([val.familyID, val.probandGender, val.siblingGender, val.siblingMotherAge, val.siblingFatherAge, val.probandMotherAge, val.probandFatherAge])

    if (merge):
        if famFile == "":
            #megabaseCountMerge(path, overlap, binsize, outputPrefix)
            skip = True
        elif geneMode:
            geneCountMergeFamily(path, outputPrefix, familyData)
        else:
            megabaseCountMergeFamily(path, overlap, binsize, outputPrefix, familyData)
    else:
        #megabaseCount(path, overlap, binsize, outputPrefix)
        skip = True
    

if __name__ == '__main__':
      main(sys.argv[1:])
