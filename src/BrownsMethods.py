'''
Created on 13 Feb 2014

@author: davykavanagh
'''
import numpy as np
from scipy import stats

def calcBrownsCombinedPMasked(snpSet, genotypeArray):
    nSNPs = snpSet.nSNPs
    pValArray, adjpValArray =  snpSet.getPvalues()
    chisq = sum(-2 * np.log(pValArray))
    adjchisq = sum(-2 * np.log(adjpValArray))
    
    genotypeArray = np.ma.masked_array(genotypeArray, np.isnan(genotypeArray))
    ms = genotypeArray.mean(axis=1)[(slice(None,None,None),None)]
    datam = genotypeArray - ms
    datass = np.ma.sqrt(np.ma.sum(datam**2, axis=1))

    runningSum = 0
    for i in xrange(nSNPs-1):
        temp = np.dot(datam[i:],datam[i].T)
        d = (datass[i:]*datass[i])
        rs = temp / d
        rs = np.absolute(rs)[1:]
        runningSum += 3.25 * np.sum(rs) +  .75 * np.dot(rs, rs)
        
    sigmaSq = 4*nSNPs+2*runningSum
    E = 2*nSNPs
    df = (2*(E*E))/sigmaSq
    runningSum = sigmaSq/(2*E)
    d = chisq/runningSum
    adjd = adjchisq/runningSum
    brownsP = stats.chi2.sf(d, df)
    adjBrownsP = stats.chi2.sf(adjd, df)

    return brownsP, adjBrownsP

def calcBrownsCombinedPNonMasked(snpSet, genotypeArray):
    nSNPs = snpSet.nSNPs
    pValArray, adjpValArray =  snpSet.getPvalues()
    chisq = sum(-2 * np.log(pValArray))
    adjchisq = sum(-2 * np.log(adjpValArray))
    
    ms = genotypeArray.mean(axis=1)[(slice(None,None,None),None)]
    datam = genotypeArray - ms
    datass = np.sqrt(stats.ss(datam,axis=1))

    runningSum = 0
    for i in xrange(nSNPs-1):
        temp = np.dot(datam[i:],datam[i].T)
        d = (datass[i:]*datass[i])
        rs = temp / d
        rs = np.absolute(rs)[1:]
        runningSum += 3.25 * np.sum(rs) +  .75 * np.dot(rs, rs)

    sigmaSq = 4*nSNPs+2*runningSum
    E = 2*nSNPs
    df = (2*(E*E))/sigmaSq
    runningSum = sigmaSq/(2*E)
    d = chisq/runningSum
    adjd = adjchisq/runningSum
    brownsP = stats.chi2.sf(d, df)
    adjBrownsP = stats.chi2.sf(adjd, df)

    return brownsP, adjBrownsP


def calcBrownsCombinedPnanRemove(snpSet, genotypeArray):
    
    nSNPs = snpSet.nSNPs
    pValArray, adjpValArray =  snpSet.getPvalues()
    chisq = sum(-2 * np.log(pValArray))
    adjchisq = sum(-2 * np.log(adjpValArray))
    
    colsWithMissingData = np.where(np.isnan(genotypeArray))[1]
    genotypeArray = np.delete(genotypeArray, colWithMissingData, 1)
    
    ms = genotypeArray.mean(axis=1)[(slice(None,None,None),None)]
    datam = genotypeArray - ms
    datass = np.sqrt(stats.ss(datam,axis=1))

    runningSum = 0
    for i in xrange(nSNPs-1):
        temp = np.dot(datam[i:],datam[i].T)
        d = (datass[i:]*datass[i])
        rs = temp / d
        rs = np.absolute(rs)[1:]
        runningSum += 3.25 * np.sum(rs) +  .75 * np.dot(rs, rs)

    sigmaSq = 4*nSNPs+2*runningSum
    E = 2*nSNPs
    df = (2*(E*E))/sigmaSq
    runningSum = sigmaSq/(2*E)
    d = chisq/runningSum
    adjd = adjchisq/runningSum
    brownsP = stats.chi2.sf(d, df)
    adjBrownsP = stats.chi2.sf(adjd, df)

    return brownsP, adjBrownsP
