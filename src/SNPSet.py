'''
Created on 12 Feb 2014

@author: davykavanagh
'''
import sys
import numpy as np

class SNPSet(object):
    '''
    class to represent sets of sets of SNPs
    '''
    def __init__(self, name, snps, plinkFileObject):
              
        self.name = name
        
        snps = set(snps)
        temp = {}
        nSNPs = 0
        extraSNPs = 0
        missingSNPs = snps.intersection(plinkFileObject.missingSNPs.keys())
        snps = snps - missingSNPs
        for snp in snps:
            try:
                thisSNP = plinkFileObject.bimInfo[snp]
                temp.setdefault(thisSNP.chrm, {}).setdefault(snp, thisSNP)
                nSNPs += 1
            except KeyError:
                extraSNPs += 1  
       
        self.snps = temp
        self.nSNPs = nSNPs
        self.missingSNPs = missingSNPs
        self.nMissingSNPs = len(self.missingSNPs)
        self.extraSNPs = extraSNPs
        self.totalSNPs = self.nSNPs + self.nMissingSNPs + self.extraSNPs

        pVals = np.empty(self.nSNPs,dtype=np.float64)
        adjpVals = np.empty(self.nSNPs,dtype=np.float64)  
        minP, minPchr, minPid, minPbp, nP05, nAdjP05 = 1, None, None, None, 0, 0
        
        idx = 0
        for i, chrm in enumerate(self.snps.values()):
            for j, snp in enumerate(chrm.values()):
                pVals[idx] = snp.pValue
                adjpVals[idx] = snp.adjPvalue
                idx+=1
                
                if snp.pValue < minP:
                    minP, minPchr, minPid, minPbp = snp.pValue, snp.chrm, snp.id, snp.bp 
                if snp.pValue <= 0.05:
                    nP05 += 1
                if snp.adjPvalue <= 0.05:
                    nAdjP05 += 1
                     
        self.bestP, self.bestPchr, self.bestPid, self.bestPbp = minP, minPchr, minPid, minPbp
        self.pValues, self.adjPvalues, self.nP05, self.nAdjP05 = pVals, adjpVals, nP05, nAdjP05
        
    def __len__(self):
        return self.nSNPs
    
    def getPvalues(self):
        return self.pValues, self.adjPvalues
    
    def overidenSNPs(self, n):
        self.nSNPs = n
    
    def __str__(self):
        return 'Set {0} with {1} snps ({2} total; {3} missing; {4} extra)'.format(self.name, 
                                                                         self.nSNPs, 
                                                                         self.totalSNPs, 
                                                                         self.nMissingSNPs,
                                                                         self.extraSNPs)
            
    def __unicode__(self):
        return self.__str__()
    
    def getSNPs(self):
        allSNPs = []
        for chr in self.snps.values():
            for snp in chr.keys():
                allSNPs.append(snp)         
        return allSNPs