'''
Created on 12 Feb 2014

@author: davykavanagh
'''
import sys

class resultsFile(object):
    '''
    classdocs
    '''
    def __init__(self, fileName):
        '''
        Constructor
        '''
        self.filename = fileName
        #Check if file name was provided
        try: 
            results = open(self.filename, 'r')
        except IOError:
            print "File {0} could not be read".format(self.filename)
            sys.exit()
        except TypeError:
            print "No results file provided".format(self.filename)
            sys.exit()
        
        header = results.next().rstrip().rsplit()
        try:
            snpColNum = header.index('SNP')
            pvalColNum = header.index('P')
        except ValueError:
            print 'Results file must have header with columns "CHR", "SNP", "BP, "P", and "ANNOT" (if running gene wide tests, or using a gene based set file). Any order of columns is allowed.'
            sys.exit()
        try:
            geneColNum = header.index('ANNOT')
        except ValueError:
            geneColNum = None
        
        snpPvals = {}
        geneAnnotations = {}
        for line in results:
            fields = line.rstrip().rsplit()
            snp, pVal = fields[snpColNum], float(fields[pvalColNum])
            snpPvals[snp] = pVal
            if geneColNum is not None:
                genes = fields[geneColNum]
                genes = genes.rsplit('|')
                geneNames = [i.rsplit('(')[0] for i in genes if i != '.']
                geneAnnotations[snp] = geneNames
                
        results.close()
        
        self.snpPvals = snpPvals 
        if geneColNum is not None:
            self.geneAnnotations = geneAnnotations
        else:
            self.geneAnnotations = None
        
        print 'Read {0} association results from file {1}'.format(len(self.snpPvals), self.filename)
        
        
    def getSNPs(self):
        return set(self.snpPvals.keys())
    
    def delete(self, snp):
        if snp in self.getSNPs():
            del self.snpPvals[snp]
            if self.geneAnnotations is not None:
                del self.geneAnnotations[snp]
            else:
                pass
        else:
            pass