import numpy as np
from scipy import stats
import sys, getopt, math, binFinder, random, csv, statsmodels.stats.multitest

class vectorAnalysis():
    def __init__(self, chrom, name, start, end, fatherAge, motherAge, varPos, count, adjCount, gender, dataset, ID):
        self.start = start
        self.end = end
        self.chrom = chrom
        self.name = name
        self.fatherAge = fatherAge
        self.motherAge = motherAge
        self.variantPosition = varPos
        self.count = count
        self.adjCount = adjCount
        self.gender = gender
        self.dataset = dataset
        self.ID = ID
        self.width = 0

def resolveMismatch(i, j, probandData, siblingData, currentChrom, minimumVariantCount):
    chromList = ["chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21", "chr22", "chr23", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chrX", "chrY"]
    probandSize = len(probandData)
    siblingSize = len(siblingData)
    done = False
    #print(probandData[i].chrom, siblingData[j].chrom, probandData[i].start, siblingData[j].start)
    if (probandData[i].chrom == siblingData[j].chrom and probandData[i].start == siblingData[j].start):
        currentChrom = probandData[i].chrom
        return i, j, currentChrom, False
    else:
        #print("mismatch detected", currentChrom, probandData[i].chrom, siblingData[j].chrom, probandData[i].start, siblingData[j].start)
        if probandData[i].chrom != siblingData[j].chrom:
            #More sibling data for the chromosome so we need to move it forward
            if probandData[i].chrom != currentChrom and siblingData[j].chrom == currentChrom:
                #print("sibling should be behind", probandData[i].chrom, siblingData[j].chrom, probandData[i].start, siblingData[j].start)
                while probandData[i].chrom != siblingData[j].chrom:
                    j += 1
                    if (j == siblingSize):
                        #print("Completed looking through sibling")
                        return i, j, currentChrom, True
                #print("sibling should be caught up", probandData[i].chrom, siblingData[j].chrom, probandData[i].start, siblingData[j].start)
                currentChrom = probandData[i].chrom
            #More proband data for the chromosome so we need to move it forward
            elif probandData[i].chrom == currentChrom and siblingData[j].chrom != currentChrom:
                #print("proband should be behind", probandData[i].chrom, siblingData[j].chrom, probandData[i].start, siblingData[j].start)
                while probandData[i].chrom != siblingData[j].chrom:
                    i += 1
                    if (i == probandSize):
                        return i, j, currentChrom, True
                #print("proband should be caught up", probandData[i].chrom, siblingData[j].chrom)
                currentChrom = probandData[i].chrom
            #Both changed to a new chromosome that isnt the same as the other.
            else:
                while chromList.index(probandData[i].chrom) < chromList.index(siblingData[j].chrom):
                    i += 1
                    if (i == probandSize):
                        return i, j, currentChrom, True
                while chromList.index(probandData[i].chrom) > chromList.index(siblingData[j].chrom):
                    j += 1
                    if (j == siblingSize):
                        return i, j, currentChrom, True

        while probandData[i].start < siblingData[j].start:
            #print("proband start behind", probandData[i].start, siblingData[j].start)
            i += 1
            if (i == probandSize):
                return i, j, currentChrom, True

        while probandData[i].start > siblingData[j].start:
            #print("sibling start behind", probandData[i].start, siblingData[j].start)
            j += 1
            if (j == siblingSize):
                return i, j, currentChrom, True

        return i, j, currentChrom, False

def meanRank(geneInfo, currentList):
    ranks = []
    weights = []
    for x in currentList:
        for i in np.arange(0, len(geneInfo)):
            if x[0] == geneInfo[i][1]:
                ranks.append(i + 1)
                weights.append(np.power(1.0 / x[1], 2))
                break
    
    return np.average(ranks, weights=weights)

def basicStats(outputPrefix):

    
    fCounts = csv.writer(open(outputPrefix + "basicStats.csv", "w"))
    fCounts.writerow(["column title here"])




def binomialCounts(probandData, siblingData, knownGenes, outputPrefix):
    countResults = []
    countPvalues = []
    knownGeneList = []
    probandID = set()
    siblingID = set()
    probandMatch = dict()
    siblingMatch = dict()
    geneFilenames = knownGenes
    currentChrom = probandData[0].chrom
    minimumVariantCount = 0
    j = 0
    i = 0

    for x in np.arange(0, len(geneFilenames)):
        knownGeneList.append(readKnownGeneList(geneFilenames[x]))

    while i < len(probandData) and j < len(siblingData):
        i, j, currentChrom, done = resolveMismatch(i, j, probandData, siblingData, currentChrom, minimumVariantCount)

        if done:
            break
        else:
            #if currentChrom[-1] in {"X", "Y"}:
                #skip = True
            #if probandData[i].count == minimumVariantCount:
                #skip = True
            if probandData[i].count == minimumVariantCount and siblingData[j].count == minimumVariantCount:
                skip = True
            elif probandData[i].count < siblingData[j].count:
                skip = True
            else:
                for y in probandData[i].ID:
                    if y in probandMatch.keys():
                        if not probandData[i].name in probandMatch[y]:
                            probandMatch[y].append(probandData[i].name)
                    else:
                        probandMatch[y] = [probandData[i].name]
                    probandID.add(y)
                for z in siblingData[j].ID:
                    if z in siblingMatch.keys():
                        if not siblingData[j].name in siblingMatch[z]:
                            siblingMatch[z].append(siblingData[j].name)
                    else:
                        siblingMatch[z] = [siblingData[j].name]
                    siblingID.add(z)

                #probandWidth = int(float(probandData[i].end)) - int(float(probandData[i].start))
                #siblingWidth = int(float(siblingData[i].end)) - int(float(siblingData[i].start))
                #siblingRatio = siblingData[j].count * 1.0 / (float(siblingData[i].end) - float(siblingData[i].start))
                table = np.array([[probandData[i].count, siblingData[i].count],[probandData[i].width - probandData[i].count, siblingData[i].width - siblingData[i].count]])

                fisherOR, fisherPvalue = stats.fisher_exact(table)
                #pvalue = stats.binom_test(probandData[i].count, probandData[i].width, siblingRatio)

                found = False
                knownGeneString = ""
                for x in np.arange(0, len(knownGeneList)):
                    if probandData[i].name in knownGeneList[x]:
                        if knownGeneString == "":
                            knownGeneString += geneFilenames[x].split("/")[-1]
                            found = True
                        else:
                            knownGeneString += "|" + geneFilenames[x].split("/")[-1]
                            found = True

                if not found:
                    knownGeneString = "No"

                countPvalues.append(fisherPvalue)
                countResults.append([probandData[i].chrom, probandData[i].name, probandData[i].start, probandData[i].end, fisherOR, fisherPvalue, probandData[i].count, siblingData[j].count, knownGeneString])

        j += 1
        i += 1
            
    if len(countPvalues) > 0:
        countBon = statsmodels.stats.multitest.multipletests(countPvalues, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
        countSidak = statsmodels.stats.multitest.multipletests(countPvalues, alpha=0.05, method='sidak', is_sorted=False, returnsorted=False)
        countHolm = statsmodels.stats.multitest.multipletests(countPvalues, alpha=0.05, method='holm', is_sorted=False, returnsorted=False)  
        countFDR = statsmodels.stats.multitest.fdrcorrection(countPvalues, alpha=0.05, method='indep', is_sorted=False)

        cCounts = csv.writer(open(outputPrefix + "peopleStats.csv", "w"))
        cCounts.writerow(["ProbandPeople", "SiblingPeople"])
        cCounts.writerow([len(probandID), len(siblingID)])

        gCounts = csv.writer(open(outputPrefix + "peopleGeneStats.csv", "w"), delimiter="\t")
        gCounts.writerow(["DataSet", "ID", "UniqueGeneCount", "Genes"])
        for x in sorted(probandMatch.keys()):
            gCounts.writerow(["proband", x, len(probandMatch[x]), probandMatch[x]])
        for x in sorted(siblingMatch.keys()):
            gCounts.writerow(["sibling", x, len(siblingMatch[x]), siblingMatch[x]])

        gCounts.writerow([len(probandID), len(siblingID)])

        fCounts = csv.writer(open(outputPrefix + "countStats.csv", "w"))
        fCounts.writerow(["Chrom", "Gene", "Start", "End", "Length", "ProbandVariantCount", "SiblingVariantCount", "OddsRatio", "Pvalue", "BonPvalue", "SidakPvalue", "HolmPvalue", "FDRPvalue", "KnownGene"])
        for row in np.arange(0, len(countResults)):
            fCounts.writerow([countResults[row][0], countResults[row][1], countResults[row][2], countResults[row][3], int(float(countResults[row][3])) - int(float(countResults[row][2])), countResults[row][6], countResults[row][7], countResults[row][4], countResults[row][5], countBon[1][row], countSidak[1][row], countHolm[1][row], countFDR[1][row], countResults[row][8]])
    else:
        print("Error: No Pvalues Generated.", file=sys.stderr)

    return countResults

def readMegaBase(filename):
    results = []
    with open(filename, mode='rt') as f:
        for line in f:
            results.append(line.strip().split(','))
    return results[1:]

def readHistogram(filename):
    results = []
    with open(filename, mode='rt') as f:
        for line in f:
            results.append(int(line.strip()))
    return results[0:-3], results[-2:]


def readCount(filename):
    results = []
    with open(filename, mode='rt') as f:
        for line in f:
            vector = line.strip().split(',')
            if (len(vector) > 1):
                for x in vector[1:]:
                    #Not sure why some entries are ''. Should check this later.
                    if(x != ''):
                        results.append(int(x))
    return results

def readVectorData(filename):
    results = []
    with open(filename, mode='rt') as f:
        header = f.readline()
        for line in f:
            s = line.strip().split("\t")

            gender = []
            dataset = []
            ID = []
            if s[4] != '[]' and s[5] != '[]':
                #print(s, s[4][1:-1], s[4][1:-1].split(","))
                fatherAge = []
                motherAge = []
                varPos = []
                starts = []
                ends = []
                width = 0
                
                count = int(float(s[7]))
                adjCount = int(float(s[8]))
                for x in s[2][1:-1].split(","):
                    if x[0] == " ":
                        starts.append(int(x[1:]))
                    else:
                        starts.append(int(x))
                for x in s[3][1:-1].split(","):
                    if x[0] == " ":
                        ends.append(int(x[1:]))
                    else:
                        ends.append(int(x))
                for x in s[4][1:-1].split(","):
                    if x.find("NA") != -1:
                        continue
                    if x[0] == " ":
                        fatherAge.append(float(x[1:]))
                    else:
                        fatherAge.append(float(x))
                for x in s[5][1:-1].split(","):
                    if x.find("NA") != -1:
                        continue
                    if x[0] == " ":
                        motherAge.append(float(x[1:]))
                    else:
                        motherAge.append(float(x))
                for x in s[6][1:-1].split(","):
                    if x.find("NA") != -1:
                        continue
                    if x[0] == " ":
                        varPos.append(float(x[1:]))
                    else:
                        varPos.append(float(x))
                for x in s[9][1:-1].split(","):
                    if x == "":
                        continue
                    if x[0] == " ":
                        gender.append(x[1:])
                    else:
                        gender.append(x)
                for x in s[10][1:-1].split(","):
                    if x == "":
                        continue
                    if x[0] == " ":
                        dataset.append(x[1:])
                    else:
                        dataset.append(x)
                for x in s[11][1:-1].split(","):
                    if x == "":
                        continue
                    if x[0] == " ":
                        ID.append(int(x[1:]))
                    else:
                        ID.append(int(x))
                    
                #print([s[0], s[1], s[2], s[3], fatherAge, motherAge])
                results.append(vectorAnalysis(s[0], s[1], starts, ends, fatherAge, motherAge, varPos, count, adjCount, gender, dataset, ID))
            else:
                #print(s)
                results.append(vectorAnalysis(s[0], s[1], [s[2]], [s[3]], [], [], [], 0, 0, gender, dataset, ID))


    #Resolve Gene Widths
    for i in results:
        if len(i.start) == 1 and len(i.end) == 1:
            i.width = i.end[0] - i.start[0]
        else:
            combined = set(range(i.start[0], i.end[0] + 1))
            for j in np.arange(1, len(i.start)):
                combined.update(range(i.start[j], i.end[j] + 1))
            i.width = len(combined)



    return results

def readKnownGeneList(filename):
    results = []
    with open(filename, mode='rt') as f:
        header = f.readline()
        for line in f:
            s = line.strip().split("\t")
            results.append(s[0])
    return results

def readKnownGeneListWeights(filename):
    results = []
    weights = []
    with open(filename, mode='rt') as f:
        #header = f.readline()
        for line in f:
            s = line.strip().split(",")
            results.append(s[1])
            weights.append(np.power(1.0 / int(s[6]), 2))

    return results, weights
            

def ageVectorStats(probandAgeVector, siblingAgeVector, outputPrefix):
    fatherResults = []
    motherResults = []
    fatherPvalues = []
    motherPvalues = []
    fatherTestStats = []
    motherTestStats = []
    currentChrom = probandAgeVector[0].chrom
    j = 0
    i = 0
    minimumVariantCount = 5

    while (i < len(probandAgeVector)) and (j < len(siblingAgeVector)):
        #print("premismatch check", i, j, currentChrom)
        i, j, currentChrom, done = resolveMismatch(i, j, probandAgeVector, siblingAgeVector, currentChrom, minimumVariantCount)
        #print("postmismatch check", i, j, currentChrom)
        if done:
            break
        else:
            #print(probandAgeVector[i].fatherAge, siblingAgeVector[j].fatherAge)
            if probandAgeVector[i].count < minimumVariantCount or siblingAgeVector[j].count < minimumVariantCount:
                skip = True
                #fatherResults.append([probandAgeVector[i].chrom, probandAgeVector[i].name, probandAgeVector[i].start, probandAgeVector[i].end, 99999, 1.0])
                #fatherPvalues.append(1.0)
            else:
                testStat, pvalue = stats.ks_2samp(probandAgeVector[i].fatherAge, siblingAgeVector[j].fatherAge)
                fatherPvalues.append(pvalue)
                fatherTestStats.append(testStat)
                fatherResults.append([probandAgeVector[i].chrom, probandAgeVector[i].name, probandAgeVector[i].start, probandAgeVector[i].end, testStat, pvalue, probandAgeVector[i].count, siblingAgeVector[j].count])
                
            if probandAgeVector[i].count < minimumVariantCount or siblingAgeVector[j].count < minimumVariantCount:
                skip = True
                #motherResults.append([probandAgeVector[i].chrom, probandAgeVector[i].name, probandAgeVector[i].start, probandAgeVector[i].end, 99999, 1.0])
                #motherPvalues.append(1.0)
            else:
                testStat, pvalue = stats.ks_2samp(probandAgeVector[i].motherAge, siblingAgeVector[j].motherAge)
                motherPvalues.append(pvalue)
                motherTestStats.append(testStat)
                motherResults.append([probandAgeVector[i].chrom, probandAgeVector[i].name, probandAgeVector[i].start, probandAgeVector[i].end, testStat, pvalue, probandAgeVector[i].count, siblingAgeVector[j].count])
        j += 1
        i += 1

    fatherBon = statsmodels.stats.multitest.multipletests(fatherPvalues, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
    fatherSidak = statsmodels.stats.multitest.multipletests(fatherPvalues, alpha=0.05, method='sidak', is_sorted=False, returnsorted=False)
    fatherHolm = statsmodels.stats.multitest.multipletests(fatherPvalues, alpha=0.05, method='holm', is_sorted=False, returnsorted=False)  
    fatherFDR = statsmodels.stats.multitest.fdrcorrection(fatherPvalues, alpha=0.05, method='indep', is_sorted=False)

    motherBon = statsmodels.stats.multitest.multipletests(motherPvalues, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
    motherSidak = statsmodels.stats.multitest.multipletests(motherPvalues, alpha=0.05, method='sidak', is_sorted=False, returnsorted=False)
    motherHolm = statsmodels.stats.multitest.multipletests(motherPvalues, alpha=0.05, method='holm', is_sorted=False, returnsorted=False)  
    motherFDR = statsmodels.stats.multitest.fdrcorrection(motherPvalues, alpha=0.05, method='indep', is_sorted=False)

    #print("Father Chisq Test Stat test:", stats.kstest(fatherTestStats, "chi2"))
    # normParam = { "scale": np.std(fatherTestStats), "loc": np.mean(fatherTestStats) }
    # exponParam = { "scale": 1.0 / np.mean(fatherTestStats) }
    #print("Father Norm Test Stat test:", stats.kstest(fatherTestStats, "norm", args=stats.norm(**normParam)))
    #print("Father Expon Test Stat test:", stats.kstest(fatherTestStats, "expon", args=stats.expon(**exponParam)))

    #print("Mother Chisq Test Stat test:", stats.kstest(motherTestStats, "chi2"))
    #print("Mother Norm Test Stat test:", stats.kstest(motherTestStats, "norm", args=stats.norm(**normParam)))
    #print("Mother Expon Test Stat test:", stats.kstest(motherTestStats, "expon", args=stats.expon(**exponParam)))
    
    #print("Father Mann Test Stat test:", stats.mannwhitneyu(fatherTestStats, motherTestStats))

    fCounts = csv.writer(open(outputPrefix + "fatherAgeStats.csv", "w"))
    fCounts.writerow(["Chrom", "Gene", "Start", "End", "Length", "ProbandVariantCount", "SiblingVariantCount", "TestStat", "Pvalue", "BonPvalue", "SidakPvalue", "HolmPvalue", "FDRPvalue"])
    for row in np.arange(0, len(fatherResults)):
        fCounts.writerow([fatherResults[row][0], fatherResults[row][1], fatherResults[row][2], fatherResults[row][3], int(float(fatherResults[row][3])) - int(float(fatherResults[row][2])), fatherResults[row][6], fatherResults[row][7], fatherResults[row][4], fatherResults[row][5], fatherBon[1][row], fatherSidak[1][row], fatherHolm[1][row], fatherFDR[1][row]])

    mCounts = csv.writer(open(outputPrefix + "motherAgeStats.csv", "w"))
    mCounts.writerow(["Chrom", "Gene", "Start", "End", "Length", "ProbandVariantCount", "SiblingVariantCount", "TestStat", "Pvalue", "BonPvalue", "SidakPvalue", "HolmPvalue", "FDRPvalue"])
    for row in np.arange(0, len(motherResults)):
        mCounts.writerow([motherResults[row][0], motherResults[row][1], motherResults[row][2], motherResults[row][3], int(float(motherResults[row][3])) - int(float(motherResults[row][2])), motherResults[row][6], motherResults[row][7], motherResults[row][4], motherResults[row][5], motherBon[1][row], motherSidak[1][row], motherHolm[1][row], motherFDR[1][row]])

    return fatherResults, motherResults


def geneCountStats(probandData, siblingData, knownGenes, outputPrefix):
    positionResults = []
    positionPvalues = []
    positionTestStats = []
    geneFilenames = knownGenes
    knownGeneList = []
    currentChrom = probandData[0].chrom
    minimumVariantCount = 5
    i = 0
    j = 0

    for x in np.arange(0, len(geneFilenames)):
        knownGeneList.append(readKnownGeneList(geneFilenames[x]))


    while i < len(probandData) and j < len(siblingData):
        i, j, currentChrom, done = resolveMismatch(i, j, probandData, siblingData, currentChrom, minimumVariantCount)
        if done:
            break
        else:
            if probandData[i].count < minimumVariantCount or siblingData[j].count < minimumVariantCount:
                skip = True
            else:
                testStat, pvalue = stats.ks_2samp(probandData[i].variantPosition, siblingData[j].variantPosition)
                found = False
                knownGeneString = ""
                for x in np.arange(0, len(knownGeneList)):
                    if probandData[i].name in knownGeneList[x]:
                        if knownGeneString == "":
                            knownGeneString += geneFilenames[x].split("/")[-1]
                            found = True
                        else:
                            knownGeneString += "|" + geneFilenames[x].split("/")[-1]
                            found = True

                if not found:
                    knownGeneString = "No"

                positionPvalues.append(pvalue)
                positionTestStats.append(testStat)
                positionResults.append([probandData[i].chrom, probandData[i].name, probandData[i].start, probandData[i].end, testStat, pvalue, probandData[i].count, siblingData[j].count, knownGeneString])
        

        j += 1
        i += 1
            
    if len(positionPvalues) > 0:
        positionBon = statsmodels.stats.multitest.multipletests(positionPvalues, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
        positionSidak = statsmodels.stats.multitest.multipletests(positionPvalues, alpha=0.05, method='sidak', is_sorted=False, returnsorted=False)
        positionHolm = statsmodels.stats.multitest.multipletests(positionPvalues, alpha=0.05, method='holm', is_sorted=False, returnsorted=False)  
        positionFDR = statsmodels.stats.multitest.fdrcorrection(positionPvalues, alpha=0.05, method='indep', is_sorted=False)


        fCounts = csv.writer(open(outputPrefix + "positionStats.csv", "w"))
        fCounts.writerow(["Chrom", "Gene", "Start", "End", "Length", "ProbandVariantCount", "SiblingVariantCount", "TestStat", "Pvalue", "BonPvalue", "SidakPvalue", "HolmPvalue", "FDRPvalue", "KnownGene"]) 
        for row in np.arange(0, len(positionResults)):
            fCounts.writerow([positionResults[row][0], positionResults[row][1], positionResults[row][2], positionResults[row][3], int(float(positionResults[row][3])) - int(float(positionResults[row][2])), positionResults[row][6], positionResults[row][7], positionResults[row][4], positionResults[row][5], positionBon[1][row], positionSidak[1][row], positionHolm[1][row], positionFDR[1][row], positionResults[row][8]])

    return positionResults


def chromosomePositionTest(probandData, siblingData, outputPrefix):
    pvalues = [[],[]]

    # for x in np.arange(0, len(probandData)):
    #     if (probandData[x].variantPosition > 4):
    #         W, pvalue = stats.shapiro(probandData[x].variantPosition)
    #         pvalue[0].append(pvalue)

    # for x in np.arange(0, len(siblingData)):
    #     W, pvalue = stats.shapiro(siblingData[x].variantPosition)
    #     pvalue[1].append(pvalue)

    # for i in np.arange(0, len(KStestResults)):
    #     resultsBon = statsmodels.stats.multitest.multipletests(pvalues[i], alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
    #     resultsSidak = statsmodels.stats.multitest.multipletests(pvalues[i], alpha=0.05, method='sidak', is_sorted=False, returnsorted=False)
    #     resultsHolm = statsmodels.stats.multitest.multipletests(pvalues[i], alpha=0.05, method='holm', is_sorted=False, returnsorted=False)  
    #     reject, adjPvalue = statsmodels.stats.multitest.fdrcorrection(pvalues[i], alpha=0.05, method='indep', is_sorted=False)

    #     wCounts = csv.writer(open(outputPrefix + output[i] + "geneNormTest.csv", "w"))
    #     wCounts.writerow(["TestStat", "Pvalue", "BonPass", "BonCorrect", "SidakPass", "SidakCorrect", "HolmsPass", "HolmsCorrect", "FDRPass", "FDRCorrect"])
    #     for row in np.arange(0, len(KStestResults[i])):
    #         wCounts.writerow([KStestResults[i][row][0], KStestResults[i][row][1], resultsBon[0][row], resultsBon[1][row], resultsSidak[0][row], resultsSidak[1][row],resultsHolm[0][row], resultsHolm[1][row], reject[row], adjPvalue[row]])


def binStats(probandData, siblingData, outputPrefix):
    output = ["count.", "mom.", "dad."]
    probandChromDict = {}
    siblingChromDict = {}
    KStestResults = [[],[],[]]
    pvalues = [[],[],[]]
    uniqueChromList = set()

    for row in probandData:
        if row[0] not in probandChromDict:
            probandChromDict[row[0]] = [[float(row[3]) / (int(row[2]) - int(row[1]))], [row[7]], [row[8]]]
            uniqueChromList.add(row[0])
        else:
            #[Count], [FatherAge], [MotherAge]
            probandChromDict[row[0]][0].append(float(row[3]) / (int(row[2]) - int(row[1])))
            probandChromDict[row[0]][1].append(row[7])
            probandChromDict[row[0]][2].append(row[8])

    for row in siblingData:
        if row[0] not in siblingChromDict:
            siblingChromDict[row[0]] = [[float(row[3]) / (int(row[2]) - int(row[1]))], [row[7]], [row[8]]]
            uniqueChromList.add(row[0])
        else:
            siblingChromDict[row[0]][0].append(float(row[3]) / (int(row[2]) - int(row[1])))
            siblingChromDict[row[0]][1].append(row[7])
            siblingChromDict[row[0]][2].append(row[8])

    for chrom in uniqueChromList:
        # if (len(chrom) > 5):
        #     continue
        if (chrom in probandChromDict.keys()) and (chrom in siblingChromDict.keys()):
            for i in np.arange(0, len(pvalues)):
                
                testStat, pvalue = stats.ks_2samp(probandChromDict[chrom][i], siblingChromDict[chrom][i])
                # fatherTestStat, fatherPvalue = stats.ks_2samp(probandChromDict[chrom][1], siblingChromDict[chrom][1])
                # motherTestStat, motherPvalue = stats.ks_2samp(probandChromDict[chrom][2], siblingChromDict[chrom][2])
                # # print("---- KS Test using Megabase Bins in Chromosome (" + chrom + ")----")
                # # print("Test Statistic:", testStat)
                # # print("P-value", pvalue, "\n")
                pvalues[i].append(pvalue)
                # pvalues[1].append(fatherPvalue)
                # pvalues[2].append(motherPvalue)
                KStestResults[i].append([chrom, testStat, pvalue])
                # KStestResults[1].append([chrom, fatherTestStat, fatherPvalue])
                # KStestResults[2].append([chrom, motherTestStat, motherPvalue])

    for i in np.arange(0, len(KStestResults)):
        resultsBon = statsmodels.stats.multitest.multipletests(pvalues[i], alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
        resultsSidak = statsmodels.stats.multitest.multipletests(pvalues[i], alpha=0.05, method='sidak', is_sorted=False, returnsorted=False)
        resultsHolm = statsmodels.stats.multitest.multipletests(pvalues[i], alpha=0.05, method='holm', is_sorted=False, returnsorted=False)  
        reject, adjPvalue = statsmodels.stats.multitest.fdrcorrection(pvalues[i], alpha=0.05, method='indep', is_sorted=False)

        wCounts = csv.writer(open(outputPrefix + output[i] + "megaBaseKS.csv", "w"))
        wCounts.writerow(["Chrom", "TestStat", "Pvalue", "BonPass", "BonCorrect", "SidakPass", "SidakCorrect", "HolmsPass", "HolmsCorrect", "FDRPass", "FDRCorrect"])
        for row in np.arange(0, len(KStestResults[i])):
            wCounts.writerow([KStestResults[i][row][0], KStestResults[i][row][1], KStestResults[i][row][2], resultsBon[0][row], resultsBon[1][row], resultsSidak[0][row], resultsSidak[1][row],resultsHolm[0][row], resultsHolm[1][row], reject[row], adjPvalue[row]])


def binStatsGene(probandData, siblingData, outputPrefix):
    output = ["count.", "mom.", "dad."]
    probandWholeGenome = [[],[],[],[]]
    siblingWholeGenome = [[],[],[],[]]
    KStestResults = [[],[],[]]
    pvalues = [[],[],[]]

    for row in probandData:
        #[Count], [FatherAge], [MotherAge]
        probandWholeGenome[0].append(row[3])
        probandWholeGenome[1].append(row[7])
        probandWholeGenome[2].append(row[8])
        probandWholeGenome[3].append(row[9])

    for row in siblingData:
        #[Count], [FatherAge], [MotherAge]
        siblingWholeGenome[0].append(row[3])
        siblingWholeGenome[1].append(row[7])
        siblingWholeGenome[2].append(row[8])
        siblingWholeGenome[3].append(row[9])

    for i in np.arange(0, len(pvalues)):      
        testStat, pvalue = stats.ks_2samp(probandWholeGenome[i], siblingWholeGenome[i])
        pvalues[i].append(pvalue)
        KStestResults[i].append([testStat, pvalue])

    for i in np.arange(0, len(KStestResults)):
        resultsBon = statsmodels.stats.multitest.multipletests(pvalues[i], alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
        resultsSidak = statsmodels.stats.multitest.multipletests(pvalues[i], alpha=0.05, method='sidak', is_sorted=False, returnsorted=False)
        resultsHolm = statsmodels.stats.multitest.multipletests(pvalues[i], alpha=0.05, method='holm', is_sorted=False, returnsorted=False)  
        reject, adjPvalue = statsmodels.stats.multitest.fdrcorrection(pvalues[i], alpha=0.05, method='indep', is_sorted=False)

        wCounts = csv.writer(open(outputPrefix + output[i] + "genomeWide.csv", "w"))
        wCounts.writerow(["TestStat", "Pvalue", "BonPass", "BonCorrect", "SidakPass", "SidakCorrect", "HolmsPass", "HolmsCorrect", "FDRPass", "FDRCorrect"])
        for row in np.arange(0, len(KStestResults[i])):
            wCounts.writerow([KStestResults[i][row][0], KStestResults[i][row][1], resultsBon[0][row], resultsBon[1][row], resultsSidak[0][row], resultsSidak[1][row],resultsHolm[0][row], resultsHolm[1][row], reject[row], adjPvalue[row]])

    # geneNormTestStatsProband, geneNormTestPvaluesProband = stats.kstest(probandWholeGenome[0], 'norm')
    # geneNormTestStatsSibling, geneNormTestPvaluesSibling = stats.kstest(siblingWholeGenome[0], 'norm')
    # print("Proband", geneNormTestStatsProband, geneNormTestPvaluesProband)
    # print("Sibling", geneNormTestStatsSibling, geneNormTestPvaluesSibling)
    
def knownGeneComparison(geneCountData, genePositionData, filenames, weightedGeneList, unrelatedGeneList, outputPrefix):
    unrelatedGeneCountPercent = []
    unrelatedGenePosPercent = []
    sortedCount = np.array(geneCountData)
    sortedPos = np.array(genePositionData)
    #sortedDad = np.array(geneFatherData)

    sortedCount = sortedCount[sortedCount[:,5].argsort()]
    sortedPos = sortedPos[sortedPos[:,5].argsort()]
    #sortedDad = sortedDad[sortedDad[:,5].argsort()]

    countSize = len(sortedCount)
    posSize = len(sortedPos)
    #dadSize = len(sortedDad)

    wCounts = csv.writer(open(outputPrefix + "countKnownGene.csv", "w"))
    wCounts.writerow(["MeanRank", "TotalRank", "PercentageRank", "pvalue", "adjustedPvalue", "genelist", "listlength", "matchedPercentage"])

    wPos = csv.writer(open(outputPrefix + "posKnownGene.csv", "w"))
    wPos.writerow(["MeanRank", "TotalRank", "PercentageRank", "pvalue", "adjustedPvalue", "genelist", "listlength", "matchedPercentage"])

    # wDad = csv.writer(open(outputPrefix + "fatherKnownGene.csv", "w"))
    # wDad.writerow(["MeanRank", "TotalRank", "PercentageRank", "pvalue", "genelist"])
    
    for x in unrelatedGeneList:
        countMean = 0
        posMean = 0
        matchedCount = 0
        matchedPos = 0
        ranks = [[],[]]
        weights = [[], []]
        compare = readKnownGeneList(x)

        for i in np.arange(0, len(compare)):
            for j in np.arange(0, len(sortedCount)):
                if sortedCount[j][1] == compare[i]:
                    ranks[0].append(j)
                    weights[0].append(1)
                    matchedCount += 1
                    break

            for j in np.arange(0, len(sortedPos)):
                if sortedPos[j][1] == compare[i]:
                    ranks[1].append(j)
                    weights[1].append(1)
                    matchedPos += 1
                    break

        if len(ranks[0]) != 0:
            countMean = np.mean(ranks[0])
        if len(ranks[1]) != 0:
            posMean = np.mean(ranks[1])

        countPvalue = stats.binom_test(countMean, n=len(sortedCount), p=.50, alternative='less')
        posPvalue = stats.binom_test(posMean, n=len(sortedPos), p=.50, alternative='less')
        

        print("Average Rank of Known Genes for Count Data in", x, ":", countMean, "/", countSize, "=", countMean/countSize, "rank with pvalue", countPvalue)
        print("Average Rank of Known Genes for Position Data in ", x, ":", posMean, "/", posSize, "=", posMean/posSize, "rank with pvalue", posPvalue)

        unrelatedGeneCountPercent.append(countMean / countSize)
        unrelatedGenePosPercent.append(posMean / posSize)

        wCounts.writerow([countMean, countSize, countMean/countSize, countPvalue, "NA", x.split("/")[-1], len(compare), matchedCount / len(compare)])
        wPos.writerow([posMean, posSize, posMean/posSize, posPvalue, "NA", x.split("/")[-1], len(compare), matchedPos / len(compare)])

    bestUnrelatedGeneCountPercent = min(unrelatedGeneCountPercent)
    bestUnrelatedGenePosPercent = min(unrelatedGenePosPercent)

    for x in filenames:
        countMean = 0
        posMean = 0
        matchedCount = 0
        matchedPos = 0
        ranks = [[],[]]
        compare = readKnownGeneList(x)

        for i in np.arange(0, len(compare)):
            for j in np.arange(0, len(sortedCount)):
                if sortedCount[j][1] == compare[i]:
                    ranks[0].append(j)
                    matchedCount += 1
                    break

            for j in np.arange(0, len(sortedPos)):
                if sortedPos[j][1] == compare[i]:
                    ranks[1].append(j)
                    matchedPos += 1
                    break
                    

        if len(ranks[0]) != 0:
            countMean = np.mean(ranks[0])
        if len(ranks[1]) != 0:
            posMean = np.mean(ranks[1])

        countPvalue = stats.binom_test(countMean, n=len(sortedCount), p=.50, alternative='less')
        posPvalue = stats.binom_test(posMean, n=len(sortedPos), p=.50, alternative='less')
        adjCountPvalue = stats.binom_test(countMean, n=len(sortedCount), p=bestUnrelatedGeneCountPercent, alternative='less')
        adjPosPvalue = stats.binom_test(posMean, n=len(sortedPos), p=bestUnrelatedGenePosPercent, alternative='less')

        print("Average Rank of Known Genes for Count Data in", x, ":", countMean, "/", countSize, "=", countMean/countSize, "rank with pvalue", countPvalue)
        print("Average Rank of Known Genes for Position Data in ", x, ":", posMean, "/", posSize, "=", posMean/posSize, "rank with pvalue", posPvalue)


        wCounts.writerow([countMean, countSize, countMean/countSize, countPvalue, adjCountPvalue, x.split("/")[-1], len(compare), matchedCount / len(compare)])
        wPos.writerow([posMean, posSize, posMean/posSize, posPvalue, adjPosPvalue, x.split("/")[-1], len(compare), matchedPos / len(compare)])

    for x in weightedGeneList:
        countMean = 0
        posMean = 0
        matchedCount = 0
        matchedPos = 0
        ranks = [[],[]]
        weights = [[], []]
        compare, geneWeights = readKnownGeneListWeights(x)

        for i in np.arange(0, len(compare)):
            for j in np.arange(0, len(sortedCount)):
                if sortedCount[j][1] == compare[i]:
                    ranks[0].append(j)
                    weights[0].append(geneWeights[i])
                    matchedCount += 1
                    break

            for j in np.arange(0, len(sortedPos)):
                if sortedPos[j][1] == compare[i]:
                    ranks[1].append(j)
                    weights[1].append(geneWeights[i])
                    matchedPos += 1
                    break

        if len(ranks[0]) != 0:
            countMean = np.average(ranks[0], None, weights[0])
        if len(ranks[1]) != 0:
            posMean = np.average(ranks[1], None, weights[1])

        countPvalue = stats.binom_test(countMean, n=len(sortedCount), p=.50, alternative='less')
        posPvalue = stats.binom_test(posMean, n=len(sortedPos), p=.50, alternative='less')
        adjCountPvalue = stats.binom_test(countMean, n=len(sortedCount), p=bestUnrelatedGeneCountPercent, alternative='less')
        adjPosPvalue = stats.binom_test(posMean, n=len(sortedPos), p=bestUnrelatedGenePosPercent, alternative='less')

        #if(x.split("/")[-1] == "gene_score_all.list"):
        #    listLengthTest(sortedCount, compare, geneWeights, countMean, outputPrefix)

        print("Average Rank of Known Genes for Count Data in", x, ":", countMean, "/", countSize, "=", countMean/countSize, "rank with pvalue", countPvalue)
        print("Average Rank of Known Genes for Position Data in ", x, ":", posMean, "/", posSize, "=", posMean/posSize, "rank with pvalue", posPvalue)

        wCounts.writerow([countMean, countSize, countMean/countSize, countPvalue, adjCountPvalue, "[W]" + x.split("/")[-1], len(compare), matchedCount / len(compare)])
        wPos.writerow([posMean, posSize, posMean/posSize, posPvalue, adjPosPvalue, "[W]" + x.split("/")[-1], len(compare), matchedPos / len(compare)])



def listLengthTest(geneInfo, geneList, geneWeights, testStat, outputPrefix):
    listLengths = [50, 100, 200, 500]
    resampleResults = [[],[],[],[]]
    weightedGeneList = []
    testStats = []
    results = []
    permutationCount = 200

    for i in np.arange(0, len(geneList)):
        weightedGeneList.append([geneList[i], geneWeights[i]])

    for i in np.arange(0, len(listLengths)):
        for j in np.arange(0, permutationCount):
            currentList = random.sample(weightedGeneList, k=listLengths[i])
            currentRank = meanRank(geneInfo, currentList)
            print("Current rank for run", j, ":", currentRank)
            resampleResults[i].append(currentRank)

    for x in resampleResults:
        currentTestStat = (np.mean(x) - testStat) / (np.std(x) / np.sqrt(permutationCount))
        print("Current Test Stat at end:", currentTestStat)
        testStats.append(currentTestStat)
        pvalue = stats.norm.cdf(currentTestStat)
        if pvalue > .5:
            pvalue = 1.0 - pvalue
        results.append([np.mean(x), np.std(x), currentTestStat, 2 * pvalue])

    wCounts = csv.writer(open(outputPrefix + "listLengthTest.csv", "w"))
    wCounts.writerow(["ResampleMeanRank", "ResampleTotalRank", "ResamplePercentageRank", "sd", "lengthTestStat", "compareTestStat", "pvalue", "listlength"])

    for x in np.arange(0, len(results)):
        wCounts.writerow([results[x][0], len(geneInfo), results[x][0] / len(geneInfo), results[x][1], results[x][2], testStat, results[x][3], listLengths[x]])



def main(argv):
    opts, args = getopt.getopt(argv, "ho:", ['proHist=', 'proCount=', 'proMega=', 'proAge=', 'sibHist=', 'sibCount=', 'sibMega=', 'sibAge=', 'output=', 'knownGene=', 'weightedList=', 'unrelatedList='])
    outputPrefix = ""
    knownGeneList = []
    weightedGeneList = []
    unrelatedGeneList = []

    for opt, arg in opts:
        if opt == '--proHist':
            probandHistogramFilename = arg
        elif opt == '--sibHist':
            siblingHistogramFilename = arg
        elif opt == '--proCount':
            probandCountFilename = arg
        elif opt == '--proAge':
            probandAgeFilename = arg
        elif opt == '--sibCount':
            siblingCountFilename = arg
        elif opt == '--proMega':
            probandMegaFilename = arg
        elif opt == '--sibMega':
            siblingMegaFilename = arg
        elif opt == '--sibAge':
            siblingAgeFilename = arg
        elif opt == '--knownGene':
            knownGeneList = arg.strip().split(",")
        elif opt == '--weightedList':
            weightedGeneList = arg.strip().split(",")
        elif opt == '--unrelatedList':
            unrelatedGeneList = arg.strip().split(",")
        elif opt in ('-h'):
            print("Use", ['--proHist', '--proCount', '--sibHist', '--sibCount'], "to input filenames")
            print("-o for output file.")
            sys.exit(0)
        elif opt in ('-o', '--output'):
            outputPrefix = arg

    #probandHistogramData, probandSpecialData = readHistogram(probandHistogramFilename)
    #siblingHistogramData, siblingSpecialData = readHistogram(siblingHistogramFilename)
    
    #probandCountData = readCount(probandCountFilename)
    #siblingCountData = readCount(siblingCountFilename)

    #print(histogramData)
    #print(countsData)

    #probandCountData = [x for x in probandCountData if x not in [2000000, 1000000]]
    #siblingCountData = [x for x in siblingCountData if x not in [2000000, 1000000]]

    #probandRemovedCount = 0
    # for x in probandCountData:
    #     if (x == 2000000) or (x == 1000000):
    #         probandCountData.remove(x)
    #         probandRemovedCount += 1

    # siblingRemovedCount = 0
    # for x in siblingCountData:
    #     if (x == 2000000) or (x == 1000000):
    #         siblingCountData.remove(x)
    #         siblingRemovedCount += 1

    #print(probandCountData)
    #print(siblingCountData)
    #print(probandRemovedCount)
    #print(siblingRemovedCount)

    # pMax = np.max(probandCountData)
    # pMin = np.min(probandCountData)
    # pMean = np.mean(probandCountData)
    # pMedian = np.median(probandCountData)
    # pMode = stats.mode(probandCountData)
    # pVariance = np.var(probandCountData)
    # pSd = np.sqrt(pVariance)

    # sMax = np.max(siblingCountData)
    # sMin = np.min(siblingCountData)
    # sMean = np.mean(siblingCountData)
    # sMedian = np.median(siblingCountData)
    # sMode = stats.mode(siblingCountData)
    # sVariance = np.var(siblingCountData)
    # sSd = np.sqrt(sVariance)

    # print("Proband Min:", pMin)
    # print("Proband Max:", pMax)
    # print("Proband Mean:", pMean)
    # print("Proband Median:", pMedian)
    # print("Proband Mode:", pMode.mode)
    # print("Proband Variance:", pVariance)
    # print("Proband SD:", pSd, "\n")

    # print("Sibling Min:", sMin)
    # print("Sibling Max:", sMax)
    # print("Sibling Mean:", sMean)
    # print("Sibling Median:", sMedian)
    # print("Sibling Mode:", sMode.mode)
    # print("Sibling Variance:", sVariance)
    # print("Sibling SD:", sSd, "\n")

    # testStat, pvalue = stats.ks_2samp(probandHistogramData, siblingHistogramData)
    # print("----KS Test using Histograms----")
    # print("Test Statistic:", testStat)
    # print("P-value", pvalue, "\n")

    # testStat, pvalue = stats.ks_2samp(probandCountData, siblingCountData)
    # print("----KS Test using Counts----")
    # print("Test Statistic:", testStat)
    # print("P-value", pvalue, "\n")
    
    # #Binwise binomial test
    # print("Binwise Binomial Test")
    # for bin in np.arange(0, binFinder.binCount() - 1):
    #     probandBinCount = probandHistogramData[bin]
    #     probandN = sum(probandHistogramData)

    #     siblingBinCount = siblingHistogramData[bin]
    #     siblingN = sum(siblingHistogramData)
    #     siblingProportion = siblingBinCount * 1.0 / siblingN

    #     print(bin, probandBinCount, probandN, siblingBinCount, siblingN, siblingProportion, probandBinCount * 1.0/ probandN)
    #     binomPvalue = stats.binom_test(probandBinCount, n=probandN, p=siblingProportion)
    #     print(binFinder.binSize(bin) + ":", '%.4f' % binomPvalue, "\n")

    
    #KS Test for Megabase Paired data
    #probandMegaBaseData = readMegaBase(probandMegaFilename)
    #siblingMegaBaseData = readMegaBase(siblingMegaFilename)

    #probandVectorData = readVectorData(probandMegaFilename[0:-15] + ".geneAgeVector.csv")
    #siblingVectorData = readVectorData(siblingMegaFilename[0:-15] + ".geneAgeVector.csv")
    probandVectorData = readVectorData(probandAgeFilename)
    siblingVectorData = readVectorData(siblingAgeFilename)


    #ageVectorStats(probandVectorData, siblingVectorData, outputPrefix)
    genePositionData = geneCountStats(probandVectorData, siblingVectorData, knownGeneList + weightedGeneList, outputPrefix)
    geneCountData = binomialCounts(probandVectorData, siblingVectorData, knownGeneList + weightedGeneList, outputPrefix)

    
    # wCounts = csv.writer(open(outputPrefix + "megaBaseKS.csv", "w"))
    # wCounts.writerow(["Chrom", "TestStat", "Pvalue", "BonPass", "BonCorrect", "SidakPass", "SidakCorrect", "HolmsPass", "HolmsCorrect", "FDRPass", "FDRCorrect"])
    # for row in np.arange(0, len(KStestResults) - 1):
    #     wCounts.writerow([KStestResults[row][0], KStestResults[row][1], KStestResults[row][2], resultsBon[0][row], resultsBon[1][row], resultsSidak[0][row], resultsSidak[1][row],resultsHolm[0][row], resultsHolm[1][row], reject[row], adjPvalue[row]])



    #binStats(probandMegaBaseData, siblingMegaBaseData, outputPrefix)
    #binStatsGene(probandMegaBaseData, siblingMegaBaseData, outputPrefix)
    #geneCountStats(probandMegaBaseData, siblingMegaBaseData, outputPrefix)

    if len(knownGeneList) != 0:
        knownGeneComparison(geneCountData, genePositionData, knownGeneList, weightedGeneList, unrelatedGeneList, outputPrefix)












if __name__ == '__main__':
      main(sys.argv[1:])