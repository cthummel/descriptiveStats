import numpy as np
from scipy import stats
import sys, getopt, math, binFinder

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

    probandCountData = [x for x in probandCountData if x not in (2000000, 1000000)]
    siblingCountData = [x for x in siblingCountData if x not in (2000000, 1000000)]

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

    pMean = np.mean(probandCountData)
    pMedian = np.median(probandCountData)
    pMode = stats.mode(probandCountData)
    pVariance = np.var(probandCountData)
    pSd = np.sqrt(pVariance)

    sMean = np.mean(siblingCountData)
    sMedian = np.median(siblingCountData)
    sMode = stats.mode(siblingCountData)
    sVariance = np.var(siblingCountData)
    sSd = np.sqrt(sVariance)

    print("Proband Mean:", pMean)
    print("Proband Median:", pMedian)
    print("Proband Mode:", pMode)
    print("Proband Variance:", pVariance)
    print("Proband SD:", pSd, "\n")

    print("Sibling Mean:", sMean)
    print("Sibling Median:", sMedian)
    print("Sibling Mode:", sMode)
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
    for bin in np.arange(0, binFinder.binCount()):
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

    probandChromDict = {}
    siblingChromDict = {}
    uniqueChromList = set()

    for row in probandMegaBaseData:
        if row[0] not in probandChromDict:
            probandChromDict[row[0]] = [row[3]]
            uniqueChromList.add(row[0])
        else:
            probandChromDict[row[0]].append(row[3])

    for row in siblingMegaBaseData:
        if row[0] not in siblingChromDict:
            siblingChromDict[row[0]] = [row[3]]
            uniqueChromList.add(row[0])
        else:
            siblingChromDict[row[0]].append(row[3])

    for chrom in uniqueChromList:
        if (chrom in probandChromDict[chrom]) and (chrom in probandChromDict[chrom]):
            testStat, pvalue = stats.ks_2samp(probandChromDict[chrom], siblingChromDict[chrom])
            print("---- KS Test using Megabase Bins in Chromosome (" + chrom + ")----")
            print("Test Statistic:", testStat)
            print("P-value", pvalue, "\n")
        












if __name__ == '__main__':
      main(sys.argv[1:])