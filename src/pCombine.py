'''
Created on 12 Feb 2014

@author: davykavanagh
'''
from resultsFile import resultsFile
from plinkFiles import plinkFiles
from plinkReader import plinkReader
from SNPSet import SNPSet
from scipy import stats
import argparse
import numpy as np
import sys
import datetime
import time

def readSetsFile(filename):
    try:
        setsFile = open(filename, 'r')
    except IOError:
        print "File {0} could not be read".format(filename)
        sys.exit()
    except TypeError:
        print "No set file specified"
        sys.exit()
        
    header = setsFile.next().rstrip().rsplit()
    temp = {} 
    for line in setsFile:
        fields = line.rstrip().rsplit()
        #first column is the set name and the 2nd columns are the snp/gene name/id.
        setName, entryName, = fields[0], fields[1]
        temp.setdefault(setName, []).append(entryName)
    sets = {k:set(v) for k,v in temp.iteritems()}
    setsFile.close()
    return header, sets

def convertGeneSets2SnpSets(geneSets, plinkFileObject):
    geneToSets = {}
    for k,v in geneSets.iteritems():
        for g in v:
            geneToSets.setdefault(g, []).append(k)
    
    snpSets = {}
    snpCount = 0
    for snp in plinkFileObject.bimInfo.values():
        for g in snp.genes:
            try:                
                setsSNPisIn = geneToSets[g]
                for s in setsSNPisIn:
                    snpSets.setdefault(s, []).append(snp.id)
                snpCount += 1
            except KeyError:
                pass
    for snp, genes in plinkFileObject.missingSNPs.iteritems():
        for g in genes:
            try:                
                setsSNPisIn = geneToSets[g]
                for s in setsSNPisIn:
                    snpSets.setdefault(s, []).append(snp)
                snpCount += 1
            except KeyError:
                pass
            
    print '{0} snps in {1} sets in the converted snp sets'.format(snpCount, len(snpSets))
    return snpSets
    
    
def makeGenesSets(plinkFileObject):
    '''
    Will go through the gene info in each SNP object and then add
    that SNP to the corresponding list value, the key for which 
    is the gene name. Note that this will add a single SNP, to all
    genes that it is annotated as belonging to.
    '''
    genicSNPs = {}
    for k,v in plinkFileObject.bimInfo.iteritems():
        for g in v.genes:
            genicSNPs.setdefault(g, []).append(k)
    for k,v in plinkFileObject.missingSNPs.iteritems():
        for g in v:
            genicSNPs.setdefault(g, []).append(k)
    genicSNPs = {k:set(v) for k,v in genicSNPs.iteritems()}
    return genicSNPs

def makeListOfSNPSets(sets, genData):
    '''
    Given a dictionary of {'SET_NAME': [snp1, snp2, ..., snpN]},
    this method constructs a list of SNPSet objects
    '''
    snpSets = []
    for k,v in sets.iteritems():
        thisSet = SNPSet(k, v, genData)
        if thisSet.nSNPs > 0:
            snpSets.append(thisSet)
    return snpSets

def calcBrownsCombinedP(snpSet, genotypes):
    nSNPs = snpSet.nSNPs
    pValArray, adjpValArray =  snpSet.getPvalues()
    chisq = sum(-2 * np.log(pValArray))
    adjchisq = sum(-2 * np.log(adjpValArray))
    runningSum = 0

    for genotypeArray in genotypes:
        nSNPsOnChr = genotypeArray.shape[0]
        genotypeArray = np.ma.masked_array(genotypeArray, np.isnan(genotypeArray))
        ms = genotypeArray.mean(axis=1)[(slice(None,None,None),None)]
        datam = genotypeArray - ms
        datass = np.ma.sqrt(np.ma.sum(datam**2, axis=1))
        for i in xrange(nSNPsOnChr-1):
            temp = np.ma.dot(datam[i:],datam[i].T)
            d = (datass[i:]*datass[i])
            rs = temp / d
            rs = np.absolute(rs)[1:]
            runningSum += 3.25 * np.sum(rs) +  .75 * np.dot(rs, rs)
        #print runningSum    
    sigmaSq = 4*nSNPs+2*runningSum
    E = 2*nSNPs
    df = (2*(E*E))/sigmaSq
    runningSum = sigmaSq/(2*E)
    d = chisq/runningSum
    adjd = adjchisq/runningSum
    brownsP = stats.chi2.sf(d, df)
    adjBrownsP = stats.chi2.sf(adjd, df)

    return brownsP, adjBrownsP

def calcSimesP(snpSet):
    pValArray, adjpValArray =  snpSet.getPvalues()

    pValArray = np.sort(pValArray)
    rank = np.array(range(len(pValArray)))
    rank = rank+1
    simesP = np.min((pValArray*len(pValArray)/rank))

    adjpValArray = np.sort(adjpValArray)
    rank = np.array(range(len(adjpValArray)))
    rank = rank+1
    adjsimesP = np.min((adjpValArray*len(adjpValArray)/rank))

    return simesP, adjsimesP
    
    


def getArgs():
    parser = argparse.ArgumentParser(description = 'Combining P-values based on LD')
    parser.add_argument('--bfile', action='store', dest='plinkFileName', help='plink bed bim fam files')
    parser.add_argument('--setsfile', action='store', dest='setsFileName', help='File defining SNPsets. Can accept one of two formats: 1) Two tab delimited columns, SET_NAME, SNP_NAME. 2) Two tab delimited columns, SET_NAME, GENE_NAME. If using format 2, the association results file must also contain gene names as annotated by plink')
    parser.add_argument('--assocfile', action='store', dest='resultsFileName', help='SNP association results file. Plink results file format is the accepted format.')
    parser.add_argument('--my-lambda', action='store', dest='myLambda', help='Can provide own lambda correction.')
    parser.add_argument('--out', action='store', dest='outputFileName', help='The output file name.')
    args = parser.parse_args()
    return args

def run(args):
    print 'Reading association results file...'
    results = resultsFile(args.resultsFileName)
    
    print 'Reading plink bim and fam files...'
    genData = plinkFiles(args.plinkFileName, results)
    
    print genData.summary()
    print 'Joining Results and Genotypes...'
    genData.addPvalues(results)
    
    print genData.summary()
    
    genData.gcCorrect(args.myLambda)
    
    if args.setsFileName is None:
        '''
        If no set file is provided, then a gene wide analysis is run on the annotated results.
        If there is no ANNOT column in the results file, an error will be thrown. 
        '''
        print 'No sets file provided...\n running gene-wide test of all genes found in annotated results file.'
        genData.addGenes(results)
        sets = makeGenesSets(genData)
    else:
        header, sets = readSetsFile(args.setsFileName)
        entryCount = sum([len(v) for v in sets.values()])
        
        if 'SNP' in header:    
            print 'Detected SNP sets file...\n{0} entries in {1} sets read from the set file'.format(entryCount, len(sets))
            
        elif 'GENE' in header:
            print 'Detected Gene sets file...\n{0} genes in {1} sets read from the set file\n Converting gene sets to snp sets...'.format(entryCount, len(sets))
            genData.addGenes(results)
            sets = convertGeneSets2SnpSets(sets, genData)
                    
        else:
            print 'The sets file must have a header with the columns "SET_NAME", and either "SNP" or "GENE", with one SET to SNP/GENE pairing per line'
            sys.exit()
        

    setsToTest = makeListOfSNPSets(sets, genData)
    
    '''
    The missing SNPS dictionary has data on all SNPs in the results file,
    not found in the bim file. But we don't want to report these if they 
    weren't in any of the sets we are actually testing, so we take the 
    intersection of all snps in all sets and the missing snps and report
    only those that are present in our sets of interest. 
    '''
    missingSNPs = set(genData.missingSNPs.keys())
    allSNPsinTests = []
    for s in setsToTest:
        allSNPsinTests.extend(s.getSNPs())
    allSNPsinTests = set(allSNPsinTests)
    missingSNPs = missingSNPs.intersection(allSNPsinTests)
    
    
    try:
        outputFile = open(args.outputFileName, 'w')
    except TypeError:
        args.outputFileName = 'pCombine.output'
        outputFile = open(args.outputFileName, 'w')
    
    
    if len(missingSNPs) > 0:
        print '{0} SNPs in results file were not found in the plink bed, bim files.\n Writing this list of SNPs to the file {1}'.format(len(missingSNPs), args.outputFileName+'_missingSNPs')
        missingSNPFile = open(args.outputFileName+'_missingSNPs', 'w')
        for s in missingSNPs:
            missingSNPFile.write('{0}\n'.format(s))
    
    reader = plinkReader(genData)
    
    print 'Combining {0} p-values across {1} sets/genes.'.format(len(allSNPsinTests) - len(missingSNPs), len(setsToTest))
    header = 'Set\tnSNPs\tnSNPsMissing\tnSNPsExtra\tBestP\tBestP_ID\tbestP_CHR\tbestP_BP\tnP05\tnAdjP05\tBrownsP\tAdjBrownsP\tSimesP\tAdjSimesP\n'
    
    outputFile.write(header)
    begin = time.time()
    progress = 0
    for v in setsToTest:
        genotypes = reader.getSNPSet(v)
        brownsP, adjBrownsP = calcBrownsCombinedP(v, genotypes)
        simesP ,simesAdjP = calcSimesP(v)
        output = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n'.format(
                   v.name, v.nSNPs, v.nMissingSNPs, v.extraSNPs,
                   v.bestP, v.bestPid, v.bestPchr, v.bestPbp,
                   v.nP05, v.nAdjP05, brownsP, adjBrownsP, simesP, simesAdjP)
        outputFile.write(output)
       
        progress += 1
        sys.stdout.write('\r {0} sets/genes done; '.format(progress))
        sys.stdout.flush()
        
    finish = time.time() - begin
    outputFile.close()
    print '\nTotal time to compute: {0}\nFinished!\n\n'.format(str(datetime.timedelta(seconds = finish)))
    

if __name__=='__main__':
    args = getArgs()
    run(args)
    