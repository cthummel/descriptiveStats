import numpy as np
from scipy import stats
import sys, getopt, math


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
    opts, args = getopt.getopt(argv, "ho:", ['proHist=', 'proCount=', 'sibHist=', 'sibCount=', 'output'])
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

    probandRemovedCount = 0
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
    print("P-value", pvalue)
    








if __name__ == '__main__':
      main(sys.argv[1:])