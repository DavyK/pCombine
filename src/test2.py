'''
Created on Feb 12, 2014

@author: Davy
'''
from pCombinev2 import getArgs, run

if __name__ == '__main__':
    args = getArgs()
    
    for i in range(1000):
        x = i+1
        inFile = '../test/randomSets/randomSet_{0}.txt'.format(x)
        outFile = '../test/randomSetsResults/randomSet_{0}.txt'.format(x)
        
        args.setsFileName = inFile
        args.outputFileName = outFile
        run(args)
