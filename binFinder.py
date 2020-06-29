binSizes = ["(-inf, -2001)", "(-2000, -1001)", "(-1000, -501)", "(-500, -201)", "(-200, -101)", "(-100, -51)", 
            "(-50, -21,)", "(-20, -11)", "-10", "-9", "-8", "-7", "-6", "-5", "-4", "-3", "-2", "-1", "0", 
            "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "(11, 20)", "(21, 50)", "(51, 100)", "(101, 200)", "(201, 500)", "(501, 1000)", "(1001, 2000)", "(2001, inf)"]

def binCount():
    return len(binSizes) - 2

def findBin(x):
    #SNV bin
    if (x == 0):
        return 18

    #Insertion bins
    elif (x > 0):
        if (x == 1):
            return 19
        elif (x == 2):
            return 20
        elif (x == 3):
            return 21
        elif (x == 4):
            return 22
        elif (x == 5):
            return 23
        elif (x == 6):
            return 24
        elif (x == 7):
            return 25
        elif (x == 8):
            return 26
        elif (x == 9):
            return 27
        elif (x == 10):
            return 28
        elif (x > 10) and (x <= 20):
            return 29
        elif (x > 20) and (x <= 50):
            return 30
        elif (x > 50) and (x <= 100):
            return 31
        elif (x > 100) and (x <= 200):
            return 32
        elif (x > 200) and (x <= 500):
            return 33
        elif (x > 500) and (x <= 1000):
            return 34
        elif (x > 1000) and (x <= 2000):
            return 35
        elif (x > 2000 and x != 1000000 and x != 2000000):
            return 36

        #Mixed Data
        elif (x == 1000000):
            return 37
        #Difficult Data
        elif (x == 2000000):
            return 38

    #Deletion Bins
    elif (x < 0):
        if (x == -1):
            return 17
        elif (x == -2):
            return 16
        elif (x == -3):
            return 15
        elif (x == -4):
            return 14
        elif (x == -5):
            return 13
        elif (x == -6):
            return 12
        elif (x == -7):
            return 11
        elif (x == -8):
            return 10
        elif (x == -9):
            return 9
        elif (x == -10):
            return 8
        elif (x < -10) and (x >= -20):
            return 7
        elif (x < -20) and (x >= -50):
            return 6
        elif (x < -50) and (x >= -100):
            return 5
        elif (x < -100) and (x >= -200):
            return 4
        elif (x < -200) and (x >= -500):
            return 3
        elif (x < -500) and (x >= -1000):
            return 2
        elif (x < -1000) and (x >= -2000):
            return 1
        elif (x > -2000):
            return 0

def binSize(x):
    return binSizes[x]