'''
Created on 12 Feb 2014

@author: davykavanagh
'''
import numpy as np

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
        genotypeArray = []
        for chr in snpSet.snps.values():
            nSNPsOnChr = len(chr.values())
            onChrgenotypes = np.zeros((nSNPsOnChr, self.nPeople)) 
            for i,v in enumerate(chr.values()):
                offset = ((v.count-1)*self.nBytesPerRead)+3
                self.bed.seek(offset,0)
                self.currBytes = self.bed.read(self.nBytesPerRead)
                #print list(self.convert())
                for j, gt in enumerate(self.convert()): onChrgenotypes[i, j] = gt
            genotypeArray.append(onChrgenotypes)
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
                    raise StopIteration()