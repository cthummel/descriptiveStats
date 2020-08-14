import numpy as np
from scipy import stats
import sys, getopt, math, binFinder, csv, statsmodels.stats.multitest

class ageAnalysis():
    def __init__(self, chrom, name, start, end, fatherAge, motherAge):
        self.start = start
        self.end = end
        self.chrom = chrom
        self.name = name
        self.fatherAge = fatherAge
        self.motherAge = motherAge


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

def ageVectorStats(probandAgeVectorFile, siblingAgeVectorFile, outputPrefix):
    fatherResults = []
    motherResults = []
    probandAgeVector = []
    siblingAgeVector = []

    with open(probandAgeVectorFile, mode='rt') as f:
        header = f.readline()
        for line in f:
            s = line.strip().split("\t")
            if s[4] != '[]' and s[5] != '[]':
                #print(s, s[4][1:-1], s[4][1:-1].split(","))
                fatherAge = []
                motherAge = []
                for x in s[4][1:-1].split(","):
                    if x[0] == " ":
                        fatherAge.append(int(x[1:]))
                    elif x == "'NA'":
                        continue
                    else:
                        fatherAge.append(int(x))
                for x in s[5][1:-1].split(","):
                    if x[0] == " ":
                        motherAge.append(int(x[1:]))
                    elif x == "'NA'":
                        continue
                    else:
                        motherAge.append(int(x))
                #print([s[0], s[1], s[2], s[3], fatherAge, motherAge])
                probandAgeVector.append(ageAnalysis(s[0], s[1], s[2], s[3], fatherAge, motherAge))
            else:
                #print(s)
                probandAgeVector.append(ageAnalysis(s[0], s[1], s[2], s[3], [], []))
            

    with open(siblingAgeVectorFile, mode='rt') as f:
        header = f.readline()
        for line in f:
            s = line.strip().split("\t")
            if s[4] != '[]' and s[5] != '[]':
                #print(s, s[4][1:-1], s[4][1:-1].split(","))
                fatherAge = []
                motherAge = []
                for x in s[4][1:-1].split(","):
                    if x[0] == " ":
                        fatherAge.append(int(x[1:]))
                    elif x == "'NA'":
                        continue
                    else:
                        fatherAge.append(int(x))
                for x in s[5][1:-1].split(","):
                    if x[0] == " ":
                        motherAge.append(int(x[1:]))
                    elif x == "'NA'":
                        continue
                    else:
                        motherAge.append(int(x))
                #print([s[0], s[1], s[2], s[3], fatherAge, motherAge])
                siblingAgeVector.append(ageAnalysis(s[0], s[1], s[2], s[3], fatherAge, motherAge))
            else:
                #print(s)
                siblingAgeVector.append(ageAnalysis(s[0], s[1], s[2], s[3], [], []))

    for i in np.arange(0, len(probandAgeVector)):
        if probandAgeVector[i].fatherAge == [] or siblingAgeVector[i].fatherAge == []:
            fatherResults.append([probandAgeVector[i].chrom, probandAgeVector[i].name, probandAgeVector[i].start, probandAgeVector[i].end, 99999, 1.0])
        else:
            testStat, pvalue = stats.ks_2samp(probandAgeVector[i].fatherAge, siblingAgeVector[i].fatherAge)
            fatherResults.append([probandAgeVector[i].chrom, probandAgeVector[i].name, probandAgeVector[i].start, probandAgeVector[i].end, testStat, pvalue])
            
        if probandAgeVector[i].motherAge == [] or siblingAgeVector[i].motherAge == []:
            motherResults.append([probandAgeVector[i].chrom, probandAgeVector[i].name, probandAgeVector[i].start, probandAgeVector[i].end, 99999, 1.0])
        else:
            testStat, pvalue = stats.ks_2samp(probandAgeVector[i].motherAge, siblingAgeVector[i].motherAge)
            motherResults.append([probandAgeVector[i].chrom, probandAgeVector[i].name, probandAgeVector[i].start, probandAgeVector[i].end, testStat, pvalue])

    fCounts = csv.writer(open(outputPrefix + "fatherAgeStats.csv", "w"))
    fCounts.writerow(["Chrom", "Gene", "Start", "End", "TestStat", "Pvalue"])
    for row in np.arange(0, len(fatherResults)):
        fCounts.writerow(fatherResults[row][0], fatherResults[row][1], fatherResults[row][2], fatherResults[row][3], fatherResults[row][4], fatherResults[row][5])

    mCounts = csv.writer(open(outputPrefix + "motherAgeStats.csv", "w"))
    mCounts.writerow(["Chrom", "Gene", "Start", "End", "TestStat", "Pvalue"])
    for row in np.arange(0, len(motherResults)):
        mCounts.writerow(motherResults[row][0], motherResults[row][1], motherResults[row][2], motherResults[row][3], motherResults[row][4], motherResults[row][5])

    return fatherResults, motherResults



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
    opts, args = getopt.getopt(argv, "ho:", ['proHist=', 'proCount=', 'proMega=', 'sibHist=', 'sibCount=', 'sibMega=', 'output='])
    outputPrefix = ""

    for opt, arg in opts:
        if opt == '--proHist':
            probandHistogramFilename = arg
        elif opt == '--sibHist':
            siblingHistogramFilename = arg
        elif opt == '--proCount':
            probandCountFilename = arg
        elif opt == '--sibCount':
            siblingCountFilename = arg
        elif opt == '--proMega':
            probandMegaFilename = arg
        elif opt == '--sibMega':
            siblingMegaFilename = arg
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

    ageVectorStats("proband.merged.MM.geneAgeVector.csv", "sibling.merged.MM.geneAgeVector.csv", outputPrefix)

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
        












if __name__ == '__main__':
      main(sys.argv[1:])