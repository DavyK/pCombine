'''
Created on 12 Feb 2014

@author: davykavanagh
'''

import sys
import numpy as np
from SNP import SNP
from Person import Person
import math
from scipy import stats

class plinkFiles(object):
    '''
    classdocs
    '''
    def __init__(self, fileName, resultsObject):
        '''
        Constructor
        '''
        self.filename = fileName
        #Check if file name was provided
        if not fileName:
            print "No plink bfile stem name provided\n ... exiting!"
            sys.exit()
        
        self.bim = fileName + '.bim'
        self.fam = fileName + '.fam'
        self.bed = fileName + '.bed'

        try:
            #try to open the bim file, exit if not
            bimData = open(self.bim, 'r')
        except IOError:
            print "File {0} could not be read".format(self.bim)
            sys.exit()
            
        resultSNPs = set(resultsObject.snpPvals.keys())
            
        #parse the bim file into a dictionary of snp objects
        nSNPs = 0
        snpInfo = {}          
        for line in bimData:
            nSNPs +=1
            thisSNP = SNP(line, nSNPs)
            if thisSNP.id in resultSNPs:
                snpInfo[thisSNP.id] = thisSNP
        bimData.close()
        self.bimInfo = snpInfo
        self.nSNPs = nSNPs
        
        try:
            #try to open the fam file, exit if not
            famData = open(self.fam, 'r')
        except IOError:
            print "File {0} could not be read".format(self.fam)
            sys.exit()
        #parse the fam file into a dictionary of Person objects
        nPeople = 0
        personInfo = {}
        for line in famData:
            nPeople += 1
            thisPerson = Person(line, nPeople)
            personInfo[thisPerson.uniqueID] = thisPerson
        famData.close()
            
        self.famInfo = personInfo
        self.nPeople = nPeople
        #these variables help with reading the bed file
        self.nGenotypes = nPeople * nSNPs
        self.nAllelesPerSNP = nPeople*2
        self.bitsPerLine = int(8*math.ceil(float(self.nAllelesPerSNP)/8))
        self.nBytesPerSNP = self.bitsPerLine/8
        self.pVals = None
        self.missingSNPs = None

    def addPvalues(self, resultsObject):
        snpsInResults = resultsObject.getSNPs()
        missingSNPs = []
        for snp in snpsInResults:
            pVal = resultsObject.snpPvals[snp]
            try:
                self.bimInfo[snp].pValue = pVal
            except KeyError:
                pass

                missingSNPs.append(snp)
        missingSNPs = set(missingSNPs)
        '''
        for snp in missingSNPs:
            resultsObject.delete(snp)
        '''
        reducedBimInfo = {}
        for k, v in self.bimInfo.iteritems():
            if k in snpsInResults:
                reducedBimInfo[k] = v
                
        self.bimInfo = reducedBimInfo
        self.nSNPs = len(self.bimInfo)
        self.missingSNPs = {snp:None for snp in missingSNPs}
        
        

    def gcCorrect(self, myLambda):
        snps = self.bimInfo.keys()
        snps = [snp for snp in snps if self.bimInfo[snp].pValue]
        allPvals = [self.bimInfo[snp].pValue for snp in snps]
        
        pValArray = np.array(allPvals, dtype=np.float64)
        chi2Array = stats.chi2.isf(pValArray, 1)
        gcLambda = (np.median(chi2Array)/0.456)
        print 'Genomic Control Lambda calculated as: {0}'.format(gcLambda)
        if myLambda is None:
            print 'Correcting with calcluated lambda: {0}'.format(gcLambda)
            gcLambda = gcLambda if gcLambda > 1 else 1
            chi2Array = chi2Array/gcLambda
        else:
            myLambda = float(myLambda)
            print 'Correcting with user provided lambda: {0}'.format(myLambda)
            chi2Array = chi2Array/myLambda
        pValArray = stats.chi2.sf(chi2Array, 1)
        for i, snp in enumerate(snps):
                self.bimInfo[snp].adjPvalue = pValArray[i]
                     
    def getPvalues(self):
        pVals = np.empty(len(self.bimInfo),dtype=np.float64)
        adjpVals = np.empty(len(self.bimInfo),dtype=np.float64)
        for i, snp in enumerate(self.bimInfo.values()):
            pVals[i] = snp.pValue
            adjpVals[i] = snp.adjPvalue
        return pVals, adjpVals
    
    def addGenes(self, resultsObject):
        if resultsObject.geneAnnotations is not None:
            for snp, geneNames in resultsObject.geneAnnotations.iteritems():
                try:
                    self.bimInfo[snp].genes = geneNames
                except KeyError:
                    self.missingSNPs[snp] = geneNames
                    
        else:
            print 'Results file must have header with columns "CHR", "SNP", "BP, "P", and "ANNOT" (if running gene wide tests, or using a gene based set file). Any order of columns is allowed.'
            sys.exit()

    
    def summary(self):
        output = 'Plink Binary Ped File Summary:\n Filename:\t"{0}"\n # Individuals:\t{1}\n # SNPs:\t{2}'.format(self.filename, self.nPeople, self.nSNPs)
        return output