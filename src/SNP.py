'''
Created on 12 Feb 2014

@author: davykavanagh
'''

class SNP(object):
    '''
    class to represent simple SNP object from plink file
    '''
    def __init__(self, line, count):
        snpInfoList = line.rstrip().rsplit()
        self.chrm = snpInfoList[0]
        self.id = snpInfoList[1]
        self.bp = snpInfoList[3]
        self.count = count
        self.pValue = None
        self.adjPvalue = None
        self.genes = None
        self.missing = False
    
    def __str__(self):
        return '{0} @ {1}:{2}'.format(self.id, self.chrm, self.bp)
    
    def __unicode__(self):
        return self.__str__()
    
    def __repr__(self):
        return self.__str__()
    
    def __hash__(self):
        return(hash(str(self)))
    
    def __cmp__(self, other):
        return cmp(str(self), str(other))
