'''
Created on Apr 19, 2013

@author: davykavanagh
'''
from plinkReader import plinkReader
from plinkFiles import plinkFiles
from resultsFile import resultsFile
from pCombine import makeGenesSets, makeListOfSNPSets

if __name__=='__main__':
    results = resultsFile('../example/test.annot')
    genData = plinkFiles('../example/test', results)
    #print '{0}, {1}, {2}, {3}, {4}, {5}'.format(plinkFile.nPeople, plinkFile.nSNPs, plinkFile.nGenotypes, plinkFile.nAllelesPerSNP, plinkFile.bitsPerLine, plinkFile.nBytesPerSNP)
    genData.addPvalues(results)
    genData.addGenes(results)
    sets = makeGenesSets(genData)
    print sets
    setsToTest = makeListOfSNPSets(sets, genData)
    print setsToTest
    reader = plinkReader(genData)
    
    for snp in reader:
        print snp
        
        
    reader = plinkReader(genData)
    for v in setsToTest:
        print reader.getSNPSet(v)
