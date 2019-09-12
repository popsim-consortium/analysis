import sys
import os
import numpy as np

inDir = sys.argv[1]
print(inDir)
ratesInWins = {}
for fileName in os.listdir(inDir):
    if fileName.endswith(".txt"):
        print(fileName)
        with open(inDir + "/" + fileName, "rt") as f:
            sys.stderr.write("reading {}/{}\n".format(inDir, fileName))
            first = True
            for line in f:
                if first:
                    first = False
                else:
                    chrom, winMid, recRate = line.strip().split()[:3]
                    if chrom not in ratesInWins:
                        ratesInWins[chrom] = []
                    winMid = int(winMid)
                    recRate = float(recRate)
                    ratesInWins[chrom].append((winMid, recRate))


def getWinLenForChrom(ratesInWinsForChrom):
    prevWin = ratesInWinsForChrom[0][0]
    winLens = {}
    for win, recRates in ratesInWinsForChrom[1:]:
        winLen = win-prevWin
        if winLen in winLens:
            winLens[winLen] += 1
        else:
            winLens[winLen] = 1
        prevWin = win
    if len(winLens) != 1:
        sys.stderr.write("window lengths not consistent within chrom arms!! ARRGHHHH!\n")
    winLens = sorted(winLens.keys(), key=lambda x: winLens[x])
    return winLens[-1]


def getWinLens(ratesInWins):
    winLens = {}
    for chrom in ratesInWins:
        winLens[chrom] = getWinLenForChrom(ratesInWins[chrom])
    return winLens


winLens = getWinLens(ratesInWins)

allRates = []
for chrom in ratesInWins:
    for win, recRate in ratesInWins[chrom]:
        allRates.append(recRate)
allRates.sort()
lenCutoff = 1/np.mean(allRates) * 1e6
rateCutoff = allRates[int(len(allRates)*0.05)]
sys.stderr.write("rate cutoff: {}; length cutoff: {}\n".format(rateCutoff, lenCutoff))

for chrom in ratesInWins:
    halfWinLen = int(winLens[chrom]/2)
    mode = 0
    runLen = 0
    runStart = 1
    for winMid, recRate in ratesInWins[chrom]:
        winStart = winMid - halfWinLen
        winEnd = winMid + halfWinLen
        if mode == 1:
            if recRate <= rateCutoff:
                mode = 0
                runLen = 1
                runStart = winStart
            else:
                pass
        elif mode == 0:
            if recRate <= rateCutoff:
                runLen += 1
            else:
                if winStart-runStart >= lenCutoff:
                    print(chrom, runStart, winStart, winStart-runStart, runLen)
                mode = 1
    if mode == 0:
        if winEnd-runStart >= lenCutoff:
            print(chrom, runStart, winEnd, winEnd-runStart, runLen)
