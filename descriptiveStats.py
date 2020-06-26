import numpy as np
import scipy, sys, getopt, math


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
            for x in line.strip().split(',')[1:]:
                #Not sure why some entries are ''. Should check this later.
                if(x != ''):
                    results.append(int(x))
    return results










def main(argv):
    opts, args = getopt.getopt(argv, "h:c:o:", ['--path'])
    outputPrefix = ""

    for opt, arg in opts:
        if opt == '-h':
            histogramFilename = arg
        elif opt in ('-c', '--counts'):
            countsFilename = arg
        elif opt in ('-o', '--output'):
            outputPrefix = arg

    histogramData, specialData = readHistogram(histogramFilename)
    countsData = readCount(countsFilename)

    print(histogramData)
    print(countsData)







if __name__ == '__main__':
      main(sys.argv[1:])