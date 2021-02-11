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
    def __init__(self, name, transcriptStart, transcriptEnd, count, adjCount, insertion, deletion, snv, fatherAge, motherAge, varPos, genderList, dataset, ID):
        self.name = name
        self.transcriptStart = transcriptStart
        self.transcriptEnd = transcriptEnd
        # self.exonStart = exonStart
        # self.exonEnd = exonEnd
        self.count = count
        self.adjCount = adjCount
        self.insertion = insertion
        self.deletion = deletion
        self.snv = snv
        self.fatherAge = fatherAge
        self.motherAge = motherAge
        self.variantPosition = varPos
        self.genderList = genderList
        self.dataset = dataset
        self.ID = ID
        self.width = 0


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

def generateGeneFileSingleCategory(category):
    result = [{}, {}, {}, {}, {}, {}, {}]
    outputName = ""

    if category == "gene":
        outputName = "gene."
        with gzip.open("gencode.v34.annotation.gff3.gz", mode='rt') as f:
            for line in f:
                s = line.strip().split("\t")
                if s[0][0] == "#":
                    continue

                infoField = s[8].strip().split(";")
                
                if s[2] == "gene" and infoField[2][10:] == "protein_coding": #10: -> gene_type=
                    geneName = infoField[3][10:]
                    for i in np.arange(0, len(result)):
                        if (s[0] not in result[i].keys()):
                            result[i][s[0]] = [geneInfo(geneName, [int(s[3])], [int(s[4])], 0, 0, 0, 0, 0, [], [], [], [], [], [])]
                        else:
                            # index = -1
                            # for k in np.arange(0, len(result[i][s[0]])):
                            #     if result[i][s[0]][k].name == geneName:
                            #         index = k
                            #         break
                            # if index < 0:
                            result[i][s[0]].append(geneInfo(geneName, [int(s[3])], [int(s[4])], 0, 0, 0, 0, 0, [], [], [], [], [], []))
                            # else:
                            #     result[i][s[0]][index].transcriptStart.append(int(s[3]))
                            #     result[i][s[0]][index].transcriptEnd.append(int(s[4]))
                                

    elif category == "exon":
        outputName = "exon."
        with gzip.open("gencode.v34.annotation.gff3.gz", mode='rt') as f:
            for line in f:
                s = line.strip().split("\t")
                if s[0][0] == "#":
                    continue

                infoField = s[8].strip().split(";")

                if s[2][:3] == "exo" and infoField[4][10:] == "protein_coding":   #exon
                    geneName = infoField[5][10:]
                    for i in np.arange(0, len(result)):
                        if (s[0] not in result[i].keys()):
                            result[i][s[0]] = [geneInfo(geneName, [int(s[3])], [int(s[4])], 0, 0, 0, 0, 0, [], [], [], [], [], [])]
                        else:
                            index = -1
                            for k in np.arange(0, len(result[i][s[0]])):
                                if result[i][s[0]][k].name == geneName:
                                    index = k
                                    break
                            if index < 0:
                                result[i][s[0]].append(geneInfo(geneName, [int(s[3])], [int(s[4])], 0, 0, 0, 0, 0, [], [], [], [], [], []))
                            else:
                                result[i][s[0]][index].transcriptStart.append(int(s[3]))
                                result[i][s[0]][index].transcriptEnd.append(int(s[4]))
        
    elif category == "CDS":
        outputName = "cds."
        with gzip.open("gencode.v34.annotation.gff3.gz", mode='rt') as f:
            for line in f:
                s = line.strip().split("\t")
                if s[0][0] == "#":
                    continue

                infoField = s[8].strip().split(";")
        
                if s[2] == "CDS" and infoField[4][10:] == "protein_coding":   #cds
                    geneName = infoField[5][10:]
                    for i in np.arange(0, len(result)):
                        if (s[0] not in result[i].keys()):
                            result[i][s[0]] = [geneInfo(geneName, [int(s[3])], [int(s[4])], 0, 0, 0, 0, 0, [], [], [], [], [], [])]
                        else:
                            index = -1
                            for k in np.arange(0, len(result[i][s[0]])):
                                if result[i][s[0]][k].name == geneName:
                                    index = k
                                    break
                            if index < 0:
                                result[i][s[0]].append(geneInfo(geneName, [int(s[3])], [int(s[4])], 0, 0, 0, 0, 0, [], [], [], [], [], []))
                            else:
                                result[i][s[0]][index].transcriptStart.append(int(s[3]))
                                result[i][s[0]][index].transcriptEnd.append(int(s[4]))
    
    elif category == "transcript":
        outputName = "trans."
        with gzip.open("gencode.v34.annotation.gff3.gz", mode='rt') as f:
            for line in f:
                s = line.strip().split("\t")
                if s[0][0] == "#":
                    continue

                infoField = s[8].strip().split(";")

                if s[2][:3] == "tra" and infoField[4][10:] == "protein_coding":   #transcript
                    geneName = infoField[5][10:]
                    for i in np.arange(0, len(result)):
                        if (s[0] not in result[i].keys()):
                            result[i][s[0]] = [geneInfo(geneName, [int(s[3])], [int(s[4])], 0, 0, 0, 0, 0, [], [], [], [], [], [])]
                        else:
                            index = -1
                            for k in np.arange(0, len(result[i][s[0]])):
                                if result[i][s[0]][k].name == geneName:
                                    index = k
                                    break
                            if index < 0:
                                result[i][s[0]].append(geneInfo(geneName, [int(s[3])], [int(s[4])], 0, 0, 0, 0, 0, [], [], [], [], [], []))
                            else:
                                result[i][s[0]][index].transcriptStart.append(int(s[3]))
                                result[i][s[0]][index].transcriptEnd.append(int(s[4]))

    elif category == "3_prime_UTR":
        outputName = "3UTR."
        with gzip.open("gencode.v34.annotation.gff3.gz", mode='rt') as f:
            for line in f:
                s = line.strip().split("\t")
                if s[0][0] == "#":
                    continue

                infoField = s[8].strip().split(";")

                if s[2][:3] == "thr" and infoField[4][10:] == "protein_coding":   #three_prime_UTR
                    geneName = infoField[5][10:]
                    for i in np.arange(0, len(result)):
                        if (s[0] not in result[i].keys()):
                            result[i][s[0]] = [geneInfo(geneName, [int(s[3])], [int(s[4])], 0, 0, 0, 0, 0, [], [], [], [], [], [])]
                        else:
                            index = -1
                            for k in np.arange(0, len(result[i][s[0]])):
                                if result[i][s[0]][k].name == geneName:
                                    index = k
                                    break
                            if index < 0:
                                result[i][s[0]].append(geneInfo(geneName, [int(s[3])], [int(s[4])], 0, 0, 0, 0, 0, [], [], [], [], [], []))
                            else:
                                result[i][s[0]][index].transcriptStart.append(int(s[3]))
                                result[i][s[0]][index].transcriptEnd.append(int(s[4]))

    elif category == "5_prime_UTR":
        outputName = "5UTR."
        with gzip.open("gencode.v34.annotation.gff3.gz", mode='rt') as f:
            for line in f:
                s = line.strip().split("\t")
                if s[0][0] == "#":
                    continue

                infoField = s[8].strip().split(";")

                if s[2][:3] == "fiv" and infoField[4][10:] == "protein_coding":   #five_prime_UTR
                    geneName = infoField[5][10:]
                    for i in np.arange(0, len(result)):
                        if (s[0] not in result[i].keys()):
                            result[i][s[0]] = [geneInfo(geneName, [int(s[3])], [int(s[4])], 0, 0, 0, 0, 0, [], [], [], [], [], [])]
                        else:
                            index = -1
                            for k in np.arange(0, len(result[i][s[0]])):
                                if result[i][s[0]][k].name == geneName:
                                    index = k
                                    break
                            if index < 0:
                                result[i][s[0]].append(geneInfo(geneName, [int(s[3])], [int(s[4])], 0, 0, 0, 0, 0, [], [], [], [], [], []))
                            else:
                                result[i][s[0]][index].transcriptStart.append(int(s[3]))
                                result[i][s[0]][index].transcriptEnd.append(int(s[4]))

    elif category == "stop_codon":
        outputName = "stop."
        with gzip.open("gencode.v34.annotation.gff3.gz", mode='rt') as f:
            for line in f:
                s = line.strip().split("\t")
                if s[0][0] == "#":
                    continue

                infoField = s[8].strip().split(";")

                if s[2][:3] == "sto" and infoField[4][10:] == "protein_coding":   #stop_codon
                    geneName = infoField[5][10:]
                    for i in np.arange(0, len(result)):
                        if (s[0] not in result[i].keys()):
                            result[i][s[0]] = [geneInfo(geneName, [int(s[3])], [int(s[4])], 0, 0, 0, 0, 0, [], [], [], [], [], [])]
                        else:
                            index = -1
                            for k in np.arange(0, len(result[i][s[0]])):
                                if result[i][s[0]][k].name == geneName:
                                    index = k
                                    break
                            if index < 0:
                                result[i][s[0]].append(geneInfo(geneName, [int(s[3])], [int(s[4])], 0, 0, 0, 0, 0, [], [], [], [], [], []))
                            else:
                                result[i][s[0]][index].transcriptStart.append(int(s[3]))
                                result[i][s[0]][index].transcriptEnd.append(int(s[4]))

    #Resolve Gene Widths
    for dataset in result:
        for chrom in dataset.keys():
            for gene in dataset[chrom]:
                for i in gene.transcriptStart:
                    if len(i.start) == 1 and len(i.end) == 1:
                        i.width = i.end[0] - i.start[0]
                    else:
                        combined = set(range(i.start[0], i.end[0] + 1))
                        for j in np.arange(1, len(i.start)):
                            combined.update(range(i.start[j], i.end[j] + 1))
                        i.width = len(combined)

    return result, outputName



def generateGeneFileCategories():

    geneResult = [{}, {}, {}, {}, {}, {}, {}]
    cdsResult = [{}, {}, {}, {}, {}, {}, {}]
    exonResult = [{}, {}, {}, {}, {}, {}, {}]
    transResult = [{}, {}, {}, {}, {}, {}, {}]
    threeUTRResult = [{}, {}, {}, {}, {}, {}, {}]
    fiveUTRResult = [{}, {}, {}, {}, {}, {}, {}]
    stopResult = [{}, {}, {}, {}, {}, {}, {}]


    #generate bins
    with gzip.open("gencode.v34.annotation.gff3.gz", mode='rt') as f:
        for line in f:
            s = line.strip().split("\t")
            if s[0][0] == "#":
                continue

            infoField = s[8].strip().split(";")
             
            if s[2] == "gene" and infoField[2][10:] == "protein_coding": #10: -> gene_type=
                geneName = infoField[3][10:]
                for i in np.arange(0, len(geneResult)):
                    if (s[0] not in geneResult[i].keys()):
                        geneResult[i][s[0]] = [geneInfo(geneName, int(s[3]), int(s[4]), 0, 0, 0, 0, 0, [], [], [], [], [], [])]
                    else:
                        geneResult[i][s[0]].append(geneInfo(geneName, int(s[3]), int(s[4]), 0, 0, 0, 0, 0, [], [], [], [], [], []))

            #elif s[2][:4] in ["exo", "CDS", "tra", "thr", "fiv", "sto"] and infoField[4][10:] == "protein_coding":

            elif s[2][:4] == "exo" and infoField[4][10:] == "protein_coding":   #exon
                geneName = infoField[5][10:]
                for i in np.arange(0, len(exonResult)):
                    if (s[0] not in exonResult[i].keys()):
                        exonResult[i][s[0]] = [geneInfo(geneName, int(s[3]), int(s[4]), 0, 0, 0, 0, 0, [], [], [], [], [], [])]
                    else:
                        exonResult[i][s[0]].append(geneInfo(geneName, int(s[3]), int(s[4]), 0, 0, 0, 0, 0, [], [], [], [], [], []))

            elif s[2] == "CDS" and infoField[4][10:] == "protein_coding":       #CDS
                geneName = infoField[5][10:]
                for i in np.arange(0, len(cdsResult)):
                    if (s[0] not in cdsResult[i].keys()):
                        cdsResult[i][s[0]] = [geneInfo(geneName, int(s[3]), int(s[4]), 0, 0, 0, 0, 0, [], [], [], [], [], [])]
                    else:
                        cdsResult[i][s[0]].append(geneInfo(geneName, int(s[3]), int(s[4]), 0, 0, 0, 0, 0, [], [], [], [], [], []))

            elif s[2][:4] == "tra" and infoField[4][10:] == "protein_coding":   #transcript
                geneName = infoField[5][10:]
                for i in np.arange(0, len(transResult)):
                    if (s[0] not in transResult[i].keys()):
                        transResult[i][s[0]] = [geneInfo(geneName, int(s[3]), int(s[4]), 0, 0, 0, 0, 0, [], [], [], [], [], [])]
                    else:
                        transResult[i][s[0]].append(geneInfo(geneName, int(s[3]), int(s[4]), 0, 0, 0, 0, 0, [], [], [], [], [], []))

            elif s[2][:4] == "thr" and infoField[4][10:] == "protein_coding":   #3_prime_UTR
                geneName = infoField[5][10:]
                for i in np.arange(0, len(threeUTRResult)):
                    if (s[0] not in threeUTRResult[i].keys()):
                        threeUTRResult[i][s[0]] = [geneInfo(geneName, int(s[3]), int(s[4]), 0, 0, 0, 0, 0, [], [], [], [], [], [])]
                    else:
                        threeUTRResult[i][s[0]].append(geneInfo(geneName, int(s[3]), int(s[4]), 0, 0, 0, 0, 0, [], [], [], [], [], []))
            
            elif s[2][:4] == "fiv" and infoField[4][10:] == "protein_coding":   #five_prime_UTR
                geneName = infoField[5][10:]
                for i in np.arange(0, len(fiveUTRResult)):
                    if (s[0] not in fiveUTRResult[i].keys()):
                        fiveUTRResult[i][s[0]] = [geneInfo(geneName, int(s[3]), int(s[4]), 0, 0, 0, 0, 0, [], [], [], [], [], [])]
                    else:
                        fiveUTRResult[i][s[0]].append(geneInfo(geneName, int(s[3]), int(s[4]), 0, 0, 0, 0, 0, [], [], [], [], [], []))

            elif s[2][:4] == "sto" and infoField[4][10:] == "protein_coding":   #stop_codon
                geneName = infoField[5][10:]
                for i in np.arange(0, len(stopResult)):
                    if (s[0] not in stopResult[i].keys()):
                        stopResult[i][s[0]] = [geneInfo(geneName, int(s[3]), int(s[4]), 0, 0, 0, 0, 0, [], [], [], [], [], [])]
                    else:
                        stopResult[i][s[0]].append(geneInfo(geneName, int(s[3]), int(s[4]), 0, 0, 0, 0, 0, [], [], [], [], [], []))

    outputNames = ["gene.", "cds.", "exon.", "trans.", "3UTR.", "5UTR.", "stop."]

    return [geneResult, cdsResult, exonResult, transResult, threeUTRResult, fiveUTRResult, stopResult], outputNames 

def geneCountMergeFamily(file, outputPrefix, familyData, geneCategory):
    probandDataSet = True
    maleCount = 0
    femaleCount = 0
    fileCount = 0
    result = geneCategory
    #[MM, MF, FM, FF, male, female, full]
    probandPeopleCount = [0, 0, 0, 0, 0, 0, 0]
    siblingPeopleCount = [0, 0, 0, 0, 0, 0, 0]
    variantsPerPerson = [[], [], [], [], [], [], []]

    for root, dirs, files in os.walk(file):
        temproot = root.strip().split("/")
        #print("root", temproot)
        if temproot[len(temproot) - 2] == "s1":
            probandDataSet = False
        for filename in files:
            if(filename[-13:] == ".FINAL.vcf.gz"):
                fileCount += 1
                siblingsPaired = True
                variantDataSets = [0,0,0]
                totalVariantCount = 0
                with gzip.open(root + filename, mode='rt') as f:
                    print("Scanning variants from file:", filename)
                    currentChrom = ""
                    currentMegaBaseIndex = 0
                    currentDataSet = 0
                    currentFatherAge = 0
                    currentMotherAge = 0
                    variantMatched = False
                    gender = 4
                    full = 6
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
                                            #maleCount += 1
                                        elif (x.probandGender == "male") and (x.siblingGender == "female"):
                                            currentDataSet = 1
                                            #maleCount += 1
                                        elif (x.probandGender == "female") and (x.siblingGender == "male"):
                                            currentDataSet = 2
                                            #femaleCount += 1
                                            gender = 5
                                        elif (x.probandGender == "female") and (x.siblingGender == "female"):
                                            currentDataSet = 3
                                            #femaleCount += 1
                                            gender = 5
                                        else:
                                            #print("Extra Proband Count Watch OUT", x.probandGender, x.siblingGender)
                                            siblingsPaired = False
                                            break
                  
                                        break            
                            else:
                                for x in familyData:
                                    if x.siblingID == sample:
                                        currentFatherAge = x.siblingFatherAge
                                        currentMotherAge = x.siblingMotherAge
                                        if (x.probandGender == "male") and (x.siblingGender == "male"):
                                            currentDataSet = 0
                                            #maleCount += 1
                                        elif (x.probandGender == "male") and (x.siblingGender == "female"):
                                            currentDataSet = 1
                                            #femaleCount += 1
                                            gender = 5
                                        elif (x.probandGender == "female") and (x.siblingGender == "male"):
                                            currentDataSet = 2
                                            #maleCount += 1                                           
                                        elif (x.probandGender == "female") and (x.siblingGender == "female"):
                                            currentDataSet = 3
                                            #femaleCount += 1
                                            gender = 5
                                        else:
                                            #print("Extra Sibling Count Watch OUT", x.probandGender, x.siblingGender)
                                            siblingsPaired = False
                                            break

                                        break    
                            continue

                        if not siblingsPaired:
                            print("Skipping data for unmatched proband and sibling")
                            break

                        variantDataSets = [currentDataSet, gender, full]

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
                        totalVariantCount += variantCount

                        #Place the variant in the right megabase bin
                        #New chromosome means we add a new key to the dictionary and append a new megabase counter.
                        ##s[0] is the chromosome, s[1] is the position
                        if (len(s[0]) > 5):
                            continue
                        if (s[0] != currentChrom):
                            #currentMegaBaseEnd = result[currentDataSet][s[0]][0].transcriptEnd
                            currentChrom = s[0]

                        for gene in result[currentDataSet][currentChrom]:
                            for index in np.arange(0, len(gene.transcriptStart)):
                                if gene.transcriptStart[index] <= int(s[1]) and gene.transcriptEnd[index] >= int(s[1]):
                                    if not variantMatched:
                                        variantMatched = True
                                        if probandDataSet:
                                            if gender == 4:
                                                maleCount += 1
                                            else:
                                                femaleCount += 1
                                            probandPeopleCount[currentDataSet] += 1
                                            probandPeopleCount[gender] += 1
                                            probandPeopleCount[full] += 1
                                        else:
                                            if gender == 4:
                                                maleCount += 1
                                            else:
                                                femaleCount += 1
                                            siblingPeopleCount[currentDataSet] += 1
                                            siblingPeopleCount[gender] += 1
                                            siblingPeopleCount[full] += 1

                                    gene.count += variantCount
                                    gene.adjCount += variantCount
                                    gene.insertion += insert
                                    gene.deletion += Del
                                    gene.snv += snv
                                    gene.variantPosition.append(int(s[1]))
                                    gene.ID.append(fileCount)
                                    if currentFatherAge != "":
                                        gene.fatherAge.append(int(currentFatherAge))
                                    else:
                                        gene.fatherAge.append("NA")
                                    if currentMotherAge != "":
                                        gene.motherAge.append(int(currentMotherAge))
                                    else:
                                        gene.motherAge.append("NA")
                                    if gender == 4:
                                        gene.genderList.append("M")
                                    else:
                                        gene.genderList.append("F")
                                    if probandDataSet:
                                        gene.dataset.append("P")
                                    else:
                                        gene.dataset.append("S")
                                    break

                            if gene.transcriptStart[0] > int(s[1]):
                                break
                        
                        #for gender combined dataset
                        for gene in result[gender][currentChrom]:
                            for index in np.arange(0, len(gene.transcriptStart)):
                                if gene.transcriptStart[index] <= int(s[1]) and gene.transcriptEnd[index] >= int(s[1]):
                                    if not variantMatched:
                                        variantMatched = True
                                        if probandDataSet:
                                            if gender == 4:
                                                maleCount += 1
                                            else:
                                                femaleCount += 1
                                            probandPeopleCount[currentDataSet] += 1
                                            probandPeopleCount[gender] += 1
                                            probandPeopleCount[full] += 1
                                        else:
                                            if gender == 4:
                                                maleCount += 1
                                            else:
                                                femaleCount += 1
                                            siblingPeopleCount[currentDataSet] += 1
                                            siblingPeopleCount[gender] += 1
                                            siblingPeopleCount[full] += 1

                                    gene.count += variantCount
                                    gene.insertion += insert
                                    gene.deletion += Del
                                    gene.snv += snv
                                    gene.variantPosition.append(int(s[1]))
                                    gene.ID.append(fileCount)
                                    if currentFatherAge != "":
                                        gene.fatherAge.append(int(currentFatherAge))
                                    else:
                                        gene.fatherAge.append("NA")
                                    if currentMotherAge != "":
                                        gene.motherAge.append(int(currentMotherAge))
                                    else:
                                        gene.motherAge.append("NA")
                                    if gender == 4:
                                        gene.genderList.append("M")
                                    else:
                                        gene.genderList.append("F")
                                    if probandDataSet:
                                        gene.dataset.append("P")
                                    else:
                                        gene.dataset.append("S")
                                    break
                            if gene.transcriptStart[0] > int(s[1]):
                                break

                        #for full dataset
                        for gene in result[full][currentChrom]:
                            for index in np.arange(0, len(gene.transcriptStart)):
                                if gene.transcriptStart[index] <= int(s[1]) and gene.transcriptEnd[index] >= int(s[1]):
                                    if not variantMatched:
                                        variantMatched = True
                                        if probandDataSet:
                                            if gender == 4:
                                                maleCount += 1
                                            else:
                                                femaleCount += 1
                                            probandPeopleCount[currentDataSet] += 1
                                            probandPeopleCount[gender] += 1
                                            probandPeopleCount[full] += 1
                                        else:
                                            if gender == 4:
                                                maleCount += 1
                                            else:
                                                femaleCount += 1
                                            siblingPeopleCount[currentDataSet] += 1
                                            siblingPeopleCount[gender] += 1
                                            siblingPeopleCount[full] += 1
                                            
                                    gene.count += variantCount
                                    gene.adjCount += variantCount
                                    gene.insertion += insert
                                    gene.deletion += Del
                                    gene.snv += snv
                                    gene.variantPosition.append(int(s[1]))
                                    gene.ID.append(fileCount)
                                    if currentFatherAge != "":
                                        gene.fatherAge.append(int(currentFatherAge))
                                    else:
                                        gene.fatherAge.append("NA")
                                    if currentMotherAge != "":
                                        gene.motherAge.append(int(currentMotherAge))
                                    else:
                                        gene.motherAge.append("NA")
                                    if gender == 4:
                                        gene.genderList.append("M")
                                    else:
                                        gene.genderList.append("F")
                                    if probandDataSet:
                                        gene.dataset.append("P")
                                    else:
                                        gene.dataset.append("S")
                                    break
                            if gene.transcriptStart[0] > int(s[1]):
                                break
                            
                        #updateChromInfoDict(result[currentDataSet], currentMegaBaseIndex, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, False)
                for dataset in np.arange(0, len(variantsPerPerson)):
                    if dataset in variantDataSets:
                        variantsPerPerson[dataset].append(totalVariantCount)
                    else:
                        variantsPerPerson[dataset].append("NA")


    outputIndex = 0
    typePrefix = ["MM", "MF", "FM", "FF", "male", "female", "full"]

    # maleCounts = csv.writer(open(outputPrefix + "male" + ".geneCounts.csv", "w"))
    # maleAgeCounts = csv.writer(open(outputPrefix + "male" + ".geneAgeVector.csv", "w"), delimiter="\t")
    # femaleCounts = csv.writer(open(outputPrefix + "female" + ".geneCounts.csv", "w"))
    # femaleAgeCounts = csv.writer(open(outputPrefix + "female" + ".geneAgeVector.csv", "w"), delimiter="\t")

    # maleCounts.writerow(["Chrom", "Start", "End", "Count", "Insertions", "Deletions", "SNV", "MeanFatherAge", "MeanMotherAge", "Gene"])
    # maleAgeCounts.writerow(["Chrom", "Gene", "Start", "End", "FatherAge", "MotherAge", "VariantPosition"])
    # femaleCounts.writerow(["Chrom", "Start", "End", "Count", "Insertions", "Deletions", "SNV", "MeanFatherAge", "MeanMotherAge", "Gene"])
    # femaleAgeCounts.writerow(["Chrom", "Gene", "Start", "End", "FatherAge", "MotherAge", "VariantPosition"])

    #Fix lopsided sample count
    if maleCount > femaleCount:
        for key in result[5].keys():
            for val in result[5][key]:
                val.adjCount = val.count * 1.0 * maleCount / femaleCount
    elif femaleCount > maleCount:
        for key in result[4].keys():
            for val in result[4][key]:
                val.adjCount = val.count * 1.0 * femaleCount / maleCount


    bCounts = csv.writer(open(outputPrefix + "geneBasicData.csv", "w"), delimiter="\t")
    bCounts.writerow(["DataSet", "MMCount", "MFCount", "FMCount", "FFCount", "MaleCount", "FemaleCount", "ratio", "FullCount", "numberOfFiles"])
    bCounts.writerow(["Proband", probandPeopleCount[0], probandPeopleCount[1], probandPeopleCount[2], probandPeopleCount[3], probandPeopleCount[4], probandPeopleCount[5], maleCount / femaleCount, probandPeopleCount[6], fileCount])
    bCounts.writerow(["Sibling", siblingPeopleCount[0], siblingPeopleCount[1], siblingPeopleCount[2], siblingPeopleCount[3], siblingPeopleCount[4], siblingPeopleCount[5], maleCount / femaleCount, siblingPeopleCount[6], fileCount])

    vCounts = csv.writer(open(outputPrefix + "variantData.csv", "w"), delimiter="\t")
    vCounts.writerow(["MMCount", "MFCount", "FMCount", "FFCount", "MaleCount", "FemaleCount", "FullCount"])
    if probandDataSet:
        for i in np.arange(0, probandPeopleCount[6]):
            vCounts.writerow([variantsPerPerson[0][i], variantsPerPerson[1][i], variantsPerPerson[2][i], variantsPerPerson[3][i], variantsPerPerson[4][i], variantsPerPerson[5][i], variantsPerPerson[6][i]])
    else:
        for i in np.arange(0, siblingPeopleCount[6]):
            vCounts.writerow([variantsPerPerson[0][i], variantsPerPerson[1][i], variantsPerPerson[2][i], variantsPerPerson[3][i], variantsPerPerson[4][i], variantsPerPerson[5][i], variantsPerPerson[6][i]])

    for x in result:
        wCounts = csv.writer(open(outputPrefix + typePrefix[outputIndex] + ".geneCounts.csv", "w"))
        fAgeCounts = csv.writer(open(outputPrefix + typePrefix[outputIndex] + ".geneAgeVector.csv", "w"), delimiter="\t")
        fAgeCounts.writerow(["Chrom", "Gene", "Start", "End", "FatherAge", "MotherAge", "VariantPosition", "Count", "AdjustCount", "Gender", "DataSet", "ID"])
        wCounts.writerow(["Chrom", "Start", "End", "Count", "Insertions", "Deletions", "SNV", "MeanFatherAge", "MeanMotherAge", "Gene", "Gender", "DataSet", "ID"])
        for key in sorted(x.keys()):
            for val in x[key]:
                fAgeCounts.writerow([key, val.name, val.transcriptStart, val.transcriptEnd, val.fatherAge, val.motherAge, val.variantPosition, val.count, val.adjCount, val.genderList, val.dataset, val.ID])
                # if outputIndex < 2:
                #     maleAgeCounts.writerow([key, val.name, val.transcriptStart, val.transcriptEnd, val.fatherAge, val.motherAge, val.variantPosition])
                # else:
                #     femaleAgeCounts.writerow([key, val.name, val.transcriptStart, val.transcriptEnd, val.fatherAge, val.motherAge, val.variantPosition])

                if val.fatherAge == [] or val.motherAge == []:
                    wCounts.writerow([key, val.transcriptStart, val.transcriptEnd, val.count, val.insertion, val.deletion, val.snv, "NA", "NA", val.name, val.genderList, val.dataset, val.ID])
                    # if outputIndex < 2:
                    #     maleCounts.writerow([key, val.transcriptStart, val.transcriptEnd, val.count, val.insertion, val.deletion, val.snv, "NA", "NA", val.name])
                    # else:
                    #     femaleCounts.writerow([key, val.transcriptStart, val.transcriptEnd, val.count, val.insertion, val.deletion, val.snv, "NA", "NA", val.name])
                else:
                    wCounts.writerow([key, val.transcriptStart, val.transcriptEnd, val.count, val.insertion, val.deletion, val.snv, ageMean(val.fatherAge), ageMean(val.motherAge), val.name, val.genderList, val.dataset, val.ID])
                    # if outputIndex < 2:
                    #     maleCounts.writerow([key, val.transcriptStart, val.transcriptEnd, val.count, val.insertion, val.deletion, val.snv, ageMean(val.fatherAge), ageMean(val.motherAge), val.name])
                    # else:
                    #     femaleCounts.writerow([key, val.transcriptStart, val.transcriptEnd, val.count, val.insertion, val.deletion, val.snv, ageMean(val.fatherAge), ageMean(val.motherAge), val.name])

        outputIndex += 1

    return result

    
def megabaseCountMergeFamily(file, overlap, binsize, outputPrefix, familyData):
    megabaseSize = binsize
    maleCount = 0
    femaleCount = 0
    probandDataSet = True
    chromInfoDict = {}
    chromInfoDictMM = {}
    chromInfoDictMF = {}
    chromInfoDictFM = {}
    chromInfoDictFF = {}
    #[MM, MF, FM, FF, male, female, full]
    result = [{}, {}, {}, {}, {}, {}, {}]

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
                    fullCurrentMegaBaseIndex = 0
                    currentDataSet = 0
                    currentFatherAge = 0
                    currentMotherAge = 0
                    gender = 4
                    full = 6
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
                                            maleCount += 1
                                        elif (x.probandGender == "male") and (x.siblingGender == "female"):
                                            currentDataSet = 1
                                            maleCount += 1
                                        elif (x.probandGender == "female") and (x.siblingGender == "male"):
                                            currentDataSet = 2
                                            femaleCount += 1
                                            gender = 5
                                        elif (x.probandGender == "female") and (x.siblingGender == "female"):
                                            currentDataSet = 3
                                            femaleCount += 1
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
                                            maleCount += 1
                                        elif (x.probandGender == "male") and (x.siblingGender == "female"):
                                            currentDataSet = 1
                                            maleCount += 1
                                        elif (x.probandGender == "female") and (x.siblingGender == "male"):
                                            currentDataSet = 2
                                            femaleCount += 1
                                            gender = 5
                                        elif (x.probandGender == "female") and (x.siblingGender == "female"):
                                            currentDataSet = 3
                                            femaleCount += 1
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
                            
                            if (s[0] not in result[full]):
                                result[full][s[0]] = [megabaseInfo(0, 0 + megabaseSize, 0, 0, 0, 0, [], [], [])]    
                                fullCurrentMegaBaseEnd = 0 + megabaseSize
                            
                            currentMegaBaseIndex = 0
                            currentMegaBaseEnd = result[gender][s[0]][0].end
                            currentChrom = s[0]

                            genderCurrentMegaBaseIndex = 0
                            genderCurrentMegaBaseEnd = result[gender][s[0]][0].end
                            genderCurrentChrom = s[0]

                            fullCurrentMegaBaseIndex = 0
                            fullCurrentMegaBaseEnd = result[full][s[0]][0].end
                            fullCurrentChrom = s[0]

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
                                genderCurrentMegaBaseEnd = result[gender][s[0]][genderCurrentMegaBaseIndex].end

                        #If we find a variant that falls in the overlap for two megabases.
                        if (int(s[1]) <= genderCurrentMegaBaseEnd and int(s[1]) > genderCurrentMegaBaseEnd - (megabaseSize * overlap)):
                            if (genderCurrentMegaBaseIndex == len(result[gender][s[0]]) - 1):
                                newMegaBaseStart = genderCurrentMegaBaseEnd - (overlap * megabaseSize)
                                appendChromInfoDict(result[gender], newMegaBaseStart, megabaseSize, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, int(s[1]))
                                #chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, variantCount, insert, Del, snv))
                                updateChromInfoDict(result[gender], genderCurrentMegaBaseIndex, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, int(s[1]), False)
                            else:
                                updateChromInfoDict(result[gender], genderCurrentMegaBaseIndex, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, int(s[1]), True)
                        else:
                            updateChromInfoDict(result[gender], genderCurrentMegaBaseIndex, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, int(s[1]), False)
                            

                        #Do it again to merge properly with whole data set   
                        while (int(s[1]) > fullCurrentMegaBaseEnd):
                            #If the megabase isnt initialized yet
                            if (fullCurrentMegaBaseIndex + 1 == len(result[full][s[0]])):
                                newMegaBaseStart = fullCurrentMegaBaseEnd - (overlap * megabaseSize)
                                #chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, 0, 0, 0, 0))
                                appendChromInfoDict(result[full], newMegaBaseStart, megabaseSize, s[0], 0, 0, 0, 0, [], [], [])
                                fullCurrentMegaBaseEnd = newMegaBaseStart + megabaseSize
                                fullCurrentMegaBaseIndex += 1
                            else:
                                fullCurrentMegaBaseIndex += 1
                                fullCurrentMegaBaseEnd = result[full][s[0]][fullCurrentMegaBaseIndex].end

                        #If we find a variant that falls in the overlap for two megabases.
                        if (int(s[1]) <= fullCurrentMegaBaseEnd and int(s[1]) > fullCurrentMegaBaseEnd - (megabaseSize * overlap)):
                            if (fullCurrentMegaBaseIndex == len(result[full][s[0]]) - 1):
                                newMegaBaseStart = fullCurrentMegaBaseEnd - (overlap * megabaseSize)
                                appendChromInfoDict(result[full], newMegaBaseStart, megabaseSize, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, int(s[1]))
                                #chromInfoDict[s[0]].append(megabaseInfo(newMegaBaseStart, newMegaBaseStart + megabaseSize, variantCount, insert, Del, snv))
                                updateChromInfoDict(result[full], fullCurrentMegaBaseIndex, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, int(s[1]), False)
                            else:
                                updateChromInfoDict(result[full], fullCurrentMegaBaseIndex, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, int(s[1]), True)
                        else:
                            updateChromInfoDict(result[full], fullCurrentMegaBaseIndex, s[0], variantCount, insert, Del, snv, currentFatherAge, currentMotherAge, int(s[1]), False)
                                
                        
                        #print(chrom, int(s[1]), chromInfoDict[chrom][currentMegaBaseIndex].end - megabaseSize, chromInfoDict[chrom][currentMegaBaseIndex].end, chromInfoDict[chrom][currentMegaBaseIndex].count)

    # maleCounts = csv.writer(open(outputPrefix + "male" + ".megaBaseCounts.csv", "w"))
    # maleAgeCounts = csv.writer(open(outputPrefix + "male" + ".megaBaseAgeVector.csv", "w"), delimiter="\t")
    # femaleCounts = csv.writer(open(outputPrefix + "female" + ".megaBaseCounts.csv", "w"))
    # femaleAgeCounts = csv.writer(open(outputPrefix + "female" + ".megaBaseAgeVector.csv", "w"), delimiter="\t")

    # maleCounts.writerow(["Chrom", "Start", "End", "Count", "Insertions", "Deletions", "SNV", "MeanFatherAge", "MeanMotherAge", "Gene"])
    # maleAgeCounts.writerow(["Chrom", "Gene", "Start", "End", "FatherAge", "MotherAge", "VariantPosition"])
    # femaleCounts.writerow(["Chrom", "Start", "End", "Count", "Insertions", "Deletions", "SNV", "MeanFatherAge", "MeanMotherAge", "Gene"])
    # femaleAgeCounts.writerow(["Chrom", "Gene", "Start", "End", "FatherAge", "MotherAge", "VariantPosition"])

    #Fix lopsided sample count
    print("MaleCount:", maleCount, "FemaleCount", femaleCount)
    if maleCount > femaleCount:
        for key in result[5].keys():
            for val in result[5][key]:
                val.count = val.count * maleCount / femaleCount
    elif femaleCount > maleCount:
        for key in result[4].keys():
            for val in result[4][key]:
                val.count = val.count * femaleCount / maleCount


    outputIndex = 0
    typePrefix = ["MM", "MF", "FM", "FF", "male", "female", "full"]
    for x in result:
        wCounts = csv.writer(open(outputPrefix + typePrefix[outputIndex] + ".megaBaseCounts.csv", "w"))
        mAgeCounts = csv.writer(open(outputPrefix + typePrefix[outputIndex] + ".megaBaseAgeVector.csv", "w"), delimiter="\t")
        mAgeCounts.writerow(["Chrom", "Gene", "Start", "End", "FatherAge", "MotherAge", "VariantPosition", "AdjustCount"])
        wCounts.writerow(["Chrom", "Start", "End", "Count", "Insertions", "Deletions", "SNV", "MeanFatherAge", "MeanMotherAge", "Gene"])
        for key in sorted(x.keys()):
            for val in x[key]:
                #print(val.fatherAge)
                mAgeCounts.writerow([key, "NA", val.start, val.end, val.fatherAge, val.motherAge, val.variantPosition, val.count])
                if val.fatherAge == [] or val.motherAge == []:
                    wCounts.writerow([key, int(val.start), int(val.end), val.count, val.insertion, val.deletion, val.snv, "NA", "NA", "NA"])
                else:
                    wCounts.writerow([key, int(val.start), int(val.end), val.count, val.insertion, val.deletion, val.snv, ageMean(val.fatherAge), ageMean(val.motherAge), "NA"])
                    
        outputIndex += 1

    return result


def main(argv):
    opts, args = getopt.getopt(argv, "h p:", ['merge=', 'output=', 'overlap=', 'binsize=', 'famFile=', 'sampleFile='])
    outputPrefix = ""
    path = '-'
    merge = False
    geneMode = False
    geneType = ""
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
            if arg.isnumeric():
                binsize = int(arg)
            else:
                geneMode = True
                geneType = arg
        elif opt == '--famFile':
            famFile = arg
        elif opt == '--sampleFile':
            sampleFile = arg
        elif opt in ('-h'):
            print("Use", ['--merge', '--output', '--overlap'], "to adjust parameters")
            sys.exit(0)

    familyData = pairSiblings(famFile, sampleFile)

    outputIndex = 0
    outputType = ["MM", "MF", "FM", "FF", "probandMale", "probandFemale", "siblingMale", "siblingFemale", "full"]
    splitFamData = [[],[],[],[],[],[],[],[],[]]
    famAge = csv.writer(open(outputPrefix + "famAge.csv", "w"))
    famAge.writerow(["famID", "probandGender", "siblingGender", "siblingMotherAge", "siblingFatherAge", "probandMotherAge", "probandFatherAge", "type"])
    for val in familyData:
        if val.probandGender == "male" and val.siblingGender == "male":
            splitFamData[0].append(val)
            splitFamData[4].append(val)
            splitFamData[6].append(val)
        elif val.probandGender == "male" and val.siblingGender == "female":
            splitFamData[1].append(val)
            splitFamData[4].append(val)
            splitFamData[7].append(val)
        elif val.probandGender == "female" and val.siblingGender == "male":
            splitFamData[2].append(val)
            splitFamData[5].append(val)
            splitFamData[6].append(val)
        elif val.probandGender == "female" and val.siblingGender == "female":
            splitFamData[3].append(val)
            splitFamData[5].append(val)
            splitFamData[7].append(val)

        splitFamData[8].append(val)

    for x in splitFamData:
        for val in x:
            famAge.writerow([val.familyID, val.probandGender, val.siblingGender, val.siblingMotherAge, val.siblingFatherAge, val.probandMotherAge, val.probandFatherAge, outputType[outputIndex]])
        outputIndex += 1

    if (merge):
        if famFile == "":
            #megabaseCountMerge(path, overlap, binsize, outputPrefix)
            skip = True
        elif geneMode:
            if geneType == "all":
                geneCategories, categoryOutputTag = generateGeneFileCategories()
                for i in np.arange(0, len(geneCategories)):
                    geneCountMergeFamily(path, outputPrefix + categoryOutputTag[i], familyData, geneCategories[i])
            else:
                geneInfo, categoryOutputTag = generateGeneFileSingleCategory(geneType)
                geneCountMergeFamily(path, outputPrefix + categoryOutputTag, familyData, geneInfo)
        else:
            megabaseCountMergeFamily(path, overlap, binsize, outputPrefix, familyData)
    else:
        #megabaseCount(path, overlap, binsize, outputPrefix)
        skip = True
    

if __name__ == '__main__':
      main(sys.argv[1:])
