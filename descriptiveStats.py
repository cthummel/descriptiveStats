import numpy as np
from scipy import stats
import sys, getopt, math, binFinder, csv, statsmodels.stats.multitest

class vectorAnalysis():
    def __init__(self, chrom, name, start, end, fatherAge, motherAge, varPos):
        self.start = start
        self.end = end
        self.chrom = chrom
        self.name = name
        self.fatherAge = fatherAge
        self.motherAge = motherAge
        self.variantPosition = varPos
        self.count = len(varPos)

def resolveMismatch(i, j, probandData, siblingData, currentChrom, minimumVariantCount):
    chromList = ["chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21", "chr22", "chr23", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chrX", "chrY"]
    probandSize = len(probandData)
    siblingSize = len(siblingData)
    done = False
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




def binomialCounts(probandData, siblingData, outputPrefix):
    countResults = []
    countPvalues = []
    currentChrom = probandData[0].chrom
    minimumVariantCount = 5
    j = 0

    for i in np.arange(0, len(probandData)):
        i, j, currentChrom, done = resolveMismatch(i, j, probandData, siblingData, currentChrom, minimumVariantCount)

        if done:
            break
        else:
            if len(probandData[i].variantPosition) < minimumVariantCount or len(siblingData[j].variantPosition) < minimumVariantCount:
                skip = True
                #fatherResults.append([probandData[i].chrom, probandData[i].name, probandData[i].start, probandData[i].end, 99999, 1.0])
                #fatherPvalues.append(1.0)
            else:
                pvalue = stats.binom_test(probandData[i].count, int(float(probandData[i].end)) - int(float(probandData[i].start)), siblingData[j].count * 1.0 / (float(siblingData[i].end) - float(siblingData[i].start)))
                countPvalues.append(pvalue)
                countResults.append([probandData[i].chrom, probandData[i].name, probandData[i].start, probandData[i].end, "NA", pvalue, len(probandData[i].variantPosition), len(siblingData[j].variantPosition)])

        j += 1
            

    countBon = statsmodels.stats.multitest.multipletests(countPvalues, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
    countSidak = statsmodels.stats.multitest.multipletests(countPvalues, alpha=0.05, method='sidak', is_sorted=False, returnsorted=False)
    countHolm = statsmodels.stats.multitest.multipletests(countPvalues, alpha=0.05, method='holm', is_sorted=False, returnsorted=False)  
    countFDR = statsmodels.stats.multitest.fdrcorrection(countPvalues, alpha=0.05, method='indep', is_sorted=False)


    fCounts = csv.writer(open(outputPrefix + "countStats.csv", "w"))
    fCounts.writerow(["Chrom", "Gene", "Start", "End", "Length", "ProbandVariantCount", "SiblingVariantCount", "TestStat", "Pvalue", "BonPvalue", "SidakPvalue", "HolmPvalue", "FDRPvalue"])
    for row in np.arange(0, len(countResults)):
        fCounts.writerow([countResults[row][0], countResults[row][1], countResults[row][2], countResults[row][3], int(float(countResults[row][3])) - int(float(countResults[row][2])), countResults[row][6], countResults[row][7], countResults[row][4], countResults[row][5], countBon[1][row], countSidak[1][row], countHolm[1][row], countFDR[1][row]])

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


def sidakProcedure(alpha, m):
    return 1 - np.power(1 - alpha, 1.0 / m)

def holmBonferroni(pvalues, alpha, m):
    k = len(pvalues)
    sortedData = np.sort(pvalues)
    for x in np.arange(0, k - 1):
        if (alpha / (m + 1 - x)) < sortedData[x]:
            return np.max(x - 1, 0)
    
    return k

def readVectorData(filename):
    results = []
    with open(filename, mode='rt') as f:
        header = f.readline()
        for line in f:
            s = line.strip().split("\t")
            if s[4] != '[]' and s[5] != '[]':
                #print(s, s[4][1:-1], s[4][1:-1].split(","))
                fatherAge = []
                motherAge = []
                varPos = []
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
                #print([s[0], s[1], s[2], s[3], fatherAge, motherAge])
                results.append(vectorAnalysis(s[0], s[1], s[2], s[3], fatherAge, motherAge, varPos))
            else:
                #print(s)
                results.append(vectorAnalysis(s[0], s[1], s[2], s[3], [], [], []))
    return results
            

def ageVectorStats(probandAgeVector, siblingAgeVector, outputPrefix):
    fatherResults = []
    motherResults = []
    fatherPvalues = []
    motherPvalues = []
    currentChrom = probandAgeVector[0].chrom
    j = 0
    minimumVariantCount = 5

    for i in np.arange(0, len(probandAgeVector)):
        i, j, currentChrom, done = resolveMismatch(i, j, probandAgeVector, siblingAgeVector, currentChrom, minimumVariantCount)

        if done:
            break
        else:
            #print(probandAgeVector[i].fatherAge, siblingAgeVector[i].fatherAge)
            if len(probandAgeVector[i].fatherAge) < minimumVariantCount or len(siblingAgeVector[i].fatherAge) < minimumVariantCount:
                skip = True
                #fatherResults.append([probandAgeVector[i].chrom, probandAgeVector[i].name, probandAgeVector[i].start, probandAgeVector[i].end, 99999, 1.0])
                #fatherPvalues.append(1.0)
            else:
                testStat, pvalue = stats.ks_2samp(probandAgeVector[i].fatherAge, siblingAgeVector[i].fatherAge)
                fatherPvalues.append(pvalue)
                fatherResults.append([probandAgeVector[i].chrom, probandAgeVector[i].name, probandAgeVector[i].start, probandAgeVector[i].end, testStat, pvalue, len(probandAgeVector[i].fatherAge), len(siblingAgeVector[i].fatherAge)])
                
            if len(probandAgeVector[i].motherAge) < minimumVariantCount or len(siblingAgeVector[i].motherAge) < minimumVariantCount:
                skip = True
                #motherResults.append([probandAgeVector[i].chrom, probandAgeVector[i].name, probandAgeVector[i].start, probandAgeVector[i].end, 99999, 1.0])
                #motherPvalues.append(1.0)
            else:
                testStat, pvalue = stats.ks_2samp(probandAgeVector[i].motherAge, siblingAgeVector[i].motherAge)
                motherPvalues.append(pvalue)
                motherResults.append([probandAgeVector[i].chrom, probandAgeVector[i].name, probandAgeVector[i].start, probandAgeVector[i].end, testStat, pvalue, len(probandAgeVector[i].motherAge), len(siblingAgeVector[i].motherAge)])

    fatherBon = statsmodels.stats.multitest.multipletests(fatherPvalues, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
    fatherSidak = statsmodels.stats.multitest.multipletests(fatherPvalues, alpha=0.05, method='sidak', is_sorted=False, returnsorted=False)
    fatherHolm = statsmodels.stats.multitest.multipletests(fatherPvalues, alpha=0.05, method='holm', is_sorted=False, returnsorted=False)  
    fatherFDR = statsmodels.stats.multitest.fdrcorrection(fatherPvalues, alpha=0.05, method='indep', is_sorted=False)

    motherBon = statsmodels.stats.multitest.multipletests(motherPvalues, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
    motherSidak = statsmodels.stats.multitest.multipletests(motherPvalues, alpha=0.05, method='sidak', is_sorted=False, returnsorted=False)
    motherHolm = statsmodels.stats.multitest.multipletests(motherPvalues, alpha=0.05, method='holm', is_sorted=False, returnsorted=False)  
    motherFDR = statsmodels.stats.multitest.fdrcorrection(motherPvalues, alpha=0.05, method='indep', is_sorted=False)


    fCounts = csv.writer(open(outputPrefix + "fatherAgeStats.csv", "w"))
    fCounts.writerow(["Chrom", "Gene", "Start", "End", "Length", "ProbandVariantCount", "SiblingVariantCount", "TestStat", "Pvalue", "BonPvalue", "SidakPvalue", "HolmPvalue", "FDRPvalue"])
    for row in np.arange(0, len(fatherResults)):
        fCounts.writerow([fatherResults[row][0], fatherResults[row][1], fatherResults[row][2], fatherResults[row][3], int(float(fatherResults[row][3])) - int(float(fatherResults[row][2])), fatherResults[row][6], fatherResults[row][7], fatherResults[row][4], fatherResults[row][5], fatherBon[1][row], fatherSidak[1][row], fatherHolm[1][row], fatherFDR[1][row]])

    mCounts = csv.writer(open(outputPrefix + "motherAgeStats.csv", "w"))
    mCounts.writerow(["Chrom", "Gene", "Start", "End", "Length", "ProbandVariantCount", "SiblingVariantCount", "TestStat", "Pvalue", "BonPvalue", "SidakPvalue", "HolmPvalue", "FDRPvalue"])
    for row in np.arange(0, len(motherResults)):
        mCounts.writerow([motherResults[row][0], motherResults[row][1], motherResults[row][2], motherResults[row][3], int(float(motherResults[row][3])) - int(float(motherResults[row][2])), motherResults[row][6], motherResults[row][7], motherResults[row][4], motherResults[row][5], motherBon[1][row], motherSidak[1][row], motherHolm[1][row], motherFDR[1][row]])

    return fatherResults, motherResults


def geneCountStats(probandData, siblingData, outputPrefix):
    positionResults = []
    positionPvalues = []
    currentChrom = probandData[0].chrom
    minimumVariantCount = 5

    j = 0

    for i in np.arange(0, len(probandData)):
        i, j, currentChrom, done = resolveMismatch(i, j, probandData, siblingData, currentChrom, minimumVariantCount)

        if done:
            break
        else:
            if len(probandData[i].variantPosition) < minimumVariantCount or len(siblingData[j].variantPosition) < minimumVariantCount:
                skip = True
            else:
                testStat, pvalue = stats.ks_2samp(probandData[i].variantPosition, siblingData[j].variantPosition)
                positionPvalues.append(pvalue)
                positionResults.append([probandData[i].chrom, probandData[i].name, probandData[i].start, probandData[i].end, testStat, pvalue, len(probandData[i].variantPosition), len(siblingData[j].variantPosition)])
        

        j += 1
            

    positionBon = statsmodels.stats.multitest.multipletests(positionPvalues, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
    positionSidak = statsmodels.stats.multitest.multipletests(positionPvalues, alpha=0.05, method='sidak', is_sorted=False, returnsorted=False)
    positionHolm = statsmodels.stats.multitest.multipletests(positionPvalues, alpha=0.05, method='holm', is_sorted=False, returnsorted=False)  
    positionFDR = statsmodels.stats.multitest.fdrcorrection(positionPvalues, alpha=0.05, method='indep', is_sorted=False)


    fCounts = csv.writer(open(outputPrefix + "positionStats.csv", "w"))
    fCounts.writerow(["Chrom", "Gene", "Start", "End", "Length", "ProbandVariantCount", "SiblingVariantCount", "TestStat", "Pvalue", "BonPvalue", "SidakPvalue", "HolmPvalue", "FDRPvalue"])
    for row in np.arange(0, len(positionResults)):
        fCounts.writerow([positionResults[row][0], positionResults[row][1], positionResults[row][2], positionResults[row][3], int(float(positionResults[row][3])) - int(float(positionResults[row][2])), positionResults[row][6], positionResults[row][7], positionResults[row][4], positionResults[row][5], positionBon[1][row], positionSidak[1][row], positionHolm[1][row], positionFDR[1][row]])

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
            probandChromDict[row[0]] = [[int(row[3]) / (int(row[2]) - int(row[1]))], [row[7]], [row[8]]]
            uniqueChromList.add(row[0])
        else:
            #[Count], [FatherAge], [MotherAge]
            probandChromDict[row[0]][0].append(int(row[3]) / (int(row[2]) - int(row[1])))
            probandChromDict[row[0]][1].append(row[7])
            probandChromDict[row[0]][2].append(row[8])

    for row in siblingData:
        if row[0] not in siblingChromDict:
            siblingChromDict[row[0]] = [[int(row[3]) / (int(row[2]) - int(row[1]))], [row[7]], [row[8]]]
            uniqueChromList.add(row[0])
        else:
            siblingChromDict[row[0]][0].append(int(row[3]) / (int(row[2]) - int(row[1])))
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
    



def main(argv):
    opts, args = getopt.getopt(argv, "ho:", ['proHist=', 'proCount=', 'proMega=', 'proAge=', 'sibHist=', 'sibCount=', 'sibMega=', 'sibAge=', 'output='])
    outputPrefix = ""

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
        elif opt in ('-h'):
            print("Use", ['--proHist', '--proCount', '--sibHist', '--sibCount'], "to input filenames")
            print("-o for output file.")
            sys.exit(0)
        elif opt in ('-o', '--output'):
            outputPrefix = arg

    probandHistogramData, probandSpecialData = readHistogram(probandHistogramFilename)
    siblingHistogramData, siblingSpecialData = readHistogram(siblingHistogramFilename)
    
    probandCountData = readCount(probandCountFilename)
    siblingCountData = readCount(siblingCountFilename)

    #print(histogramData)
    #print(countsData)

    probandCountData = [x for x in probandCountData if x not in [2000000, 1000000]]
    siblingCountData = [x for x in siblingCountData if x not in [2000000, 1000000]]

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

    pMax = np.max(probandCountData)
    pMin = np.min(probandCountData)
    pMean = np.mean(probandCountData)
    pMedian = np.median(probandCountData)
    pMode = stats.mode(probandCountData)
    pVariance = np.var(probandCountData)
    pSd = np.sqrt(pVariance)

    sMax = np.max(siblingCountData)
    sMin = np.min(siblingCountData)
    sMean = np.mean(siblingCountData)
    sMedian = np.median(siblingCountData)
    sMode = stats.mode(siblingCountData)
    sVariance = np.var(siblingCountData)
    sSd = np.sqrt(sVariance)

    print("Proband Min:", pMin)
    print("Proband Max:", pMax)
    print("Proband Mean:", pMean)
    print("Proband Median:", pMedian)
    print("Proband Mode:", pMode.mode)
    print("Proband Variance:", pVariance)
    print("Proband SD:", pSd, "\n")

    print("Sibling Min:", sMin)
    print("Sibling Max:", sMax)
    print("Sibling Mean:", sMean)
    print("Sibling Median:", sMedian)
    print("Sibling Mode:", sMode.mode)
    print("Sibling Variance:", sVariance)
    print("Sibling SD:", sSd, "\n")

    testStat, pvalue = stats.ks_2samp(probandHistogramData, siblingHistogramData)
    print("----KS Test using Histograms----")
    print("Test Statistic:", testStat)
    print("P-value", pvalue, "\n")

    testStat, pvalue = stats.ks_2samp(probandCountData, siblingCountData)
    print("----KS Test using Counts----")
    print("Test Statistic:", testStat)
    print("P-value", pvalue, "\n")
    
    #Binwise binomial test
    print("Binwise Binomial Test")
    for bin in np.arange(0, binFinder.binCount() - 1):
        probandBinCount = probandHistogramData[bin]
        probandN = sum(probandHistogramData)

        siblingBinCount = siblingHistogramData[bin]
        siblingN = sum(siblingHistogramData)
        siblingProportion = siblingBinCount * 1.0 / siblingN

        print(bin, probandBinCount, probandN, siblingBinCount, siblingN, siblingProportion, probandBinCount * 1.0/ probandN)
        binomPvalue = stats.binom_test(probandBinCount, n=probandN, p=siblingProportion)
        print(binFinder.binSize(bin) + ":", '%.4f' % binomPvalue, "\n")

    
    #KS Test for Megabase Paired data
    probandMegaBaseData = readMegaBase(probandMegaFilename)
    siblingMegaBaseData = readMegaBase(siblingMegaFilename)

    #probandVectorData = readVectorData(probandMegaFilename[0:-15] + ".geneAgeVector.csv")
    #siblingVectorData = readVectorData(siblingMegaFilename[0:-15] + ".geneAgeVector.csv")
    probandVectorData = readVectorData(probandAgeFilename)
    siblingVectorData = readVectorData(siblingAgeFilename)


    ageVectorStats(probandVectorData, siblingVectorData, outputPrefix)
    geneCountStats(probandVectorData, siblingVectorData, outputPrefix)
    binomialCounts(probandVectorData, siblingVectorData, outputPrefix)

    # probandChromDict = {}
    # siblingChromDict = {}
    # KStestResults = []
    # pvalues = []
    # uniqueChromList = set()

    # for row in probandMegaBaseData:
    #     if row[0] not in probandChromDict:
    #         probandChromDict[row[0]] = [row[3]]
    #         uniqueChromList.add(row[0])
    #     else:
    #         probandChromDict[row[0]].append(row[3])

    # for row in siblingMegaBaseData:
    #     if row[0] not in siblingChromDict:
    #         siblingChromDict[row[0]] = [row[3]]
    #         uniqueChromList.add(row[0])
    #     else:
    #         siblingChromDict[row[0]].append(row[3])

    # for chrom in uniqueChromList:
    #     # if (len(chrom) > 5):
    #     #     continue
    #     if (chrom in probandChromDict.keys()) and (chrom in siblingChromDict.keys()):
    #         testStat, pvalue = stats.ks_2samp(probandChromDict[chrom], siblingChromDict[chrom])
    #         print("---- KS Test using Megabase Bins in Chromosome (" + chrom + ")----")
    #         print("Test Statistic:", testStat)
    #         print("P-value", pvalue, "\n")
    #         pvalues.append(pvalue)
    #         KStestResults.append([chrom, testStat, pvalue])

    # # binBinomialTestResults = []
    # # for chrom in uniqueChromList:
    # #     if (chrom in probandChromDict.keys()) and (chrom in siblingChromDict.keys()):


    
    # resultsBon = statsmodels.stats.multitest.multipletests(pvalues, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)
    # resultsSidak = statsmodels.stats.multitest.multipletests(pvalues, alpha=0.05, method='sidak', is_sorted=False, returnsorted=False)
    # resultsHolm = statsmodels.stats.multitest.multipletests(pvalues, alpha=0.05, method='holm', is_sorted=False, returnsorted=False)  
    # reject, adjPvalue = statsmodels.stats.multitest.fdrcorrection(pvalues, alpha=0.05, method='indep', is_sorted=False)

    # #print(resultsBon, resultsHolm, resultsSidak)

    # wCounts = csv.writer(open(outputPrefix + "megaBaseKS.csv", "w"))
    # wCounts.writerow(["Chrom", "TestStat", "Pvalue", "BonPass", "BonCorrect", "SidakPass", "SidakCorrect", "HolmsPass", "HolmsCorrect", "FDRPass", "FDRCorrect"])
    # for row in np.arange(0, len(KStestResults) - 1):
    #     wCounts.writerow([KStestResults[row][0], KStestResults[row][1], KStestResults[row][2], resultsBon[0][row], resultsBon[1][row], resultsSidak[0][row], resultsSidak[1][row],resultsHolm[0][row], resultsHolm[1][row], reject[row], adjPvalue[row]])



    binStats(probandMegaBaseData, siblingMegaBaseData, outputPrefix)
    binStatsGene(probandMegaBaseData, siblingMegaBaseData, outputPrefix)
    #geneCountStats(probandMegaBaseData, siblingMegaBaseData, outputPrefix)


    # tempProband = [34296170, 34038118, 33969149, 33804533, 33897533, 34262966, 34101972, 34302323, 34046133, 34499237, 34264576, 34403749, 34042175, 33556121, 34418068, 34373311, 34428879, 34428880, 34004365, 33686874, 34023799, 34422787, 34249934, 33920501, 33907808, 34318097]
    # tempSibling = [33533207, 34304953, 33711308, 33554223, 33710378, 33759165, 33688167, 33931893, 34151661, 34151797, 33639300, 33639301, 33501920, 33725904, 34330898, 33707310, 33768958, 34099877, 33589316, 33548425, 33707411, 34156030]
    # print(stats.ks_2samp(tempProband, tempSibling))
    
        












if __name__ == '__main__':
      main(sys.argv[1:])