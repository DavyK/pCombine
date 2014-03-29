'''
Created on 12 Feb 2014

@author: davykavanagh
'''

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