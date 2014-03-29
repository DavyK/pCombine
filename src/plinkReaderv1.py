'''
Created on 7 Feb 2014

@author: davykavanagh
'''
import sys
import math
import numpy as np
from scipy import stats

class plinkFiles(object):
    '''
    classdocs
    '''
    def __init__(self, fileName):
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
        #parse the bim file into a dictionary of snp objects
        nSNPs = 0
        snpInfo = {}          
        for line in bimData:
            nSNPs +=1
            thisSNP = SNP(line, nSNPs)
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

    def addPvalues(self, resultsFile):
        try:
            results = open(resultsFile, 'r')
        except IOError:
            print "File {0} could not be read".format(resultsFile)
            sys.exit()
        except TypeError:
            print "No association results file specified"
            sys.exit()
        
        header = results.next().rstrip().rsplit()
        try:
            chrColNum = header.index('CHR')
            snpColNum = header.index('SNP')
            bpColNum = header.index('BP')
            pvalColNum = header.index('P')
        except ValueError:
            print 'Results file must have header with columns "CHR", "SNP", "BP and "P" (any order is allowed).'
            sys.exit()
            
        for line in results:
            fields = line.rstrip().rsplit()
            snp, pVal = fields[snpColNum], fields[pvalColNum]
            try:
                self.bimInfo[snp].pValue = float(pVal)
            except KeyError:
                chrm = fields[chrColNum]
                bp = fields[bpColNum]
                line = '{0}\t{1}\t0\t{2}\tNA\tNA\n'.format(chrm, snp, bp)
                missingSNP = SNP(line, 0)
                missingSNP.pValue = float(pVal)
                missingSNP.missing = True
                self.bimInfo[snp] = missingSNP
        results.close()

    def gcCorrect(self, myLambda):
        snps = self.bimInfo.keys()
        allPvals = [self.bimInfo[snp].pValue for snp in snps if self.bimInfo[snp].pValue]
        pValArray = np.array(allPvals, dtype=np.float64)
        chi2Array = stats.chi2.isf(pValArray, 1)
        if myLambda is None:
            gcLambda = (np.median(chi2Array)/0.456)
            print 'Genomic Control Lambda calculated as: {0}'.format(gcLambda)
            gcLambda = gcLambda if gcLambda > 1 else 1
            chi2Array = chi2Array/gcLambda
        else:
            myLambda = float(myLambda)
            print 'Correcting with user provided lambda: {0}'.format(myLambda)
            chi2Array = chi2Array/myLambda
        pValArray = stats.chi2.sf(chi2Array, 1)
        for i, snp in enumerate(snps):
                self.bimInfo[snp].adjPvalue = pValArray[i]
        sys.exit()
                     
    def getPvalues(self):
        pVals = np.empty(len(self.bimInfo),dtype=np.float64)
        adjpVals = np.empty(len(self.bimInfo),dtype=np.float64)
        for i, snp in enumerate(self.bimInfo.values()):
            pVals[i] = snp.pValue
            adjpVals[i] = snp.adjPvalue
        return pVals, adjpVals
    
    def addGenes(self, annotatedResultsFile):
        try:
            annotations = open(annotatedResultsFile, 'r')
        except IOError:
            print "File {0} could not be read".format(annotatedResultsFile)
            sys.exit()
        except TypeError:
            print "No association results file specified"
            sys.exit()
            
        header = annotations.next().rstrip().rsplit()
        try:
            snpColNum = header.index('SNP')
            geneColNum = header.index('ANNOT')
        except ValueError:
            print 'File must have header with columns ""CHR", "SNP", "BP", "P", and "ANNOT" (any order is allowed).'
            sys.exit()
        
        for line in annotations:
            fields = line.rstrip().rsplit()
            snp, genes = fields[snpColNum], fields[geneColNum]
            genes = genes.rsplit('|')
            geneNames = [i.rsplit('(')[0] for i in genes if i != '.']
            self.bimInfo[snp].genes = geneNames
        annotations.close()          
    
    def summary(self):
        output = 'Plink Binary Ped File Summary:\n Filename:\t"{0}"\n # Individuals:\t{1}\n # SNPs:\t{2}'.format(self.filename, self.nPeople, self.nSNPs)
        return output
        
class plinkReader(object):
    
    def __init__(self, plinkFilesObject):
        self.bed = open(plinkFilesObject.bed, 'rb')
        self.bed.read(3)
        self.nBytesPerRead = plinkFilesObject.nBytesPerSNP
        self.nAllelesPerSNP = plinkFilesObject.nAllelesPerSNP
        self.nPeople = plinkFilesObject.nPeople
    
    def __iter__(self):
        return self
    
    def next(self):
        self.currBytes = self.bed.read(self.nBytesPerRead)
        
        if self.currBytes == '':
            self.bed.close()
            raise StopIteration
        else:
            genotypeArray = np.zeros(self.nPeople)
            for i, gt in enumerate(self.convert()): genotypeArray[i] = gt 
            return genotypeArray
        
    def getSNPSet(self, snpSet):
        genotypeArray = np.zeros((len(snpSet), self.nPeople))
        for i, v in enumerate(snpSet.snps.values()):
            offset = ((v.count*2)-1)+3
            self.bed.seek(offset)
            self.currBytes = self.bed.read(self.nBytesPerRead)
            for j, gt in enumerate(self.convert()): genotypeArray[i, j] = gt 
        return genotypeArray
                     
    def convert(self):
        nBits = self.nAllelesPerSNP
        for byte in self.currBytes:
            for p in range(4):
                bits = (ord(byte) >> p*2) & 3
                gt = None if bits == 1 else 1 if bits == 2 else 2 if bits == 3 else 0
                yield gt
                nBits -= 2
                if nBits <= 0:
                    raise StopIteration
    

class SNPSet(object):
    '''
    class to represent sets of sets of SNPs
    '''
    def __init__(self, name, snps, plinkFileObject):
        self.name = name
        temp = {}
        snpsNotFound = []
        for snp in snps:
            thisSNP = plinkFileObject.bimInfo[snp]
            if thisSNP.missing:
                snpsNotFound.append(thisSNP)
            elif thisSNP.pValue is not None:
                temp[snp] = thisSNP
                
        self.snps = temp
        self.nSNPs = len(self.snps)
        self.missingSNPs = snpsNotFound
        self.nMissingSNPs = len(self.missingSNPs)
        self.totalSNPs = self.nSNPs + self.nMissingSNPs
        
        pVals = np.empty(self.nSNPs,dtype=np.float64)
        adjpVals = np.empty(self.nSNPs,dtype=np.float64)
        minP = 1
        minPchr = None
        minPid = None
        minPbp = None
        nP05 = 0
        minAdjP = 1
        minAdjPchr = None
        minAdjPid = None
        minAdjPbp = None
        nAdjP05 = 0
        for i, snp in enumerate(self.snps.values()):
            pVals[i] = snp.pValue
            adjpVals[i] = snp.adjPvalue 
            if snp.pValue < minP:
                minP, minPchr, minPid, minPbp = snp.pValue, snp.chr, snp.id, snp.bp
                
            if snp.adjPvalue < minAdjP:
                minAdjP, minAdjPchr, minAdjPid, minAdjPbp = snp.adjPvalue, snp.chr, snp.id, snp.bp
                
            if snp.pValue <= 0.05:
                nP05 += 1
                 
            if snp.adjPvalue <= 0.05:
                nAdjP05 += 1 
        self.bestP, self.bestPchr, self.bestPid, self.bestPbp = minP, minPchr, minPid, minPbp
        self.bestAdjP, self.bestAdjPchr, self.bestAdjPid, self.bestAdjPbp = minAdjP, minAdjPchr, minAdjPid, minAdjPbp
        self.pValues = pVals
        self.adjPvalues = adjpVals
        self.nP05 = nP05
        self.nAdjP05 = nAdjP05
        
        
    def __len__(self):
        return self.nSNPs
    
    def getPvalues(self):
        return self.pValues, self.adjPvalues
    
    def overidenSNPs(self, n):
        self.nSNPs = n
    
    def __str__(self):
        return 'Set {0} with {1} snps ({2} present; {3} missing)'.format(self.name, 
                                                                         self.totalSNPs, 
                                                                         self.nSNPs, 
                                                                         self.nMissingSNPs)
            
    def __unicode__(self):
        return self.__str__()            
            
class SNP(object):
    '''
    class to represent simple SNP object from plink file
    '''
    def __init__(self, line, count):
        snpInfoList = line.rstrip().rsplit()
        self.chr = snpInfoList[0]
        self.id = snpInfoList[1]
        self.bp = snpInfoList[3]
        self.count = count
        self.pValue = None
        self.adjPvalue = None
        self.genes = None
        self.missing = False
    
    def __str__(self):
        return '{0} @ {1}:{2} [{3}/{4}]'.format(self.id, self.chr, self.bp, self.allele1, self.allele2)
    
    def __unicode__(self):
        return self.__str__()
    
    def __hash__(self):
        return(hash(str(self)))
    
    def __cmp__(self, other):
        return cmp(str(self), str(other))

class Person(object):
    '''
        basic class to represent an individual from the fam file
    '''
    def __init__(self, line, count):
        famInfoList = line.rstrip().rsplit()
        self.FID = famInfoList[0]
        self.IID = famInfoList[1]
        self.uniqueID = self.FID + self.IID
        self.PID = famInfoList[2]
        self.MID = famInfoList[3]
        self.SEX = famInfoList[4]
        self.PHENO = famInfoList[5]
        self.count = count
        
if __name__=='__main__':
    plinkFile = plinkFiles('../test/chr22_test')
    print '{0}, {1}, {2}, {3}, {4}, {5}'.format(plinkFile.nPeople, plinkFile.nSNPs, 
                                                plinkFile.nGenotypes, plinkFile.nAllelesPerSNP, 
                                                plinkFile.bitsPerLine, plinkFile.nBytesPerSNP)