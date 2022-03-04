import pandas as pd
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
import pdb
import scipy
from scipy.optimize import minimize
from scipy.stats import multivariate_normal
from tqdm.notebook import tqdm
from scipy.io import loadmat
import math
from utils_YRange import *

def fittingStatistics_v13(responses_all, 
                         sem_all,
                         amplitude,
                        noLightSigmoidLowRSquared, 
                        midLightSigmoidLowRSquared, 
                        highLightSigmoidLowRSquared, 
                        noLightGaussLowRSquared, 
                        midLightGaussLowRSquared, 
                        highLightGaussLowRSquared,
                        noLightSigmoidLowRSquaredFitParams, 
                        midLightSigmoidLowRSquaredFitParams, 
                        highLightSigmoidLowRSquaredFitParams, 
                        noLightGaussLowRSquaredFitParams,
                        midLightGaussLowRSquaredFitParams, 
                        highLightGaussLowRSquaredFitParams,
                        noLightSigmoidLowRSquaredValues, 
                        midLightSigmoidLowRSquaredValues, 
                        highLightSigmoidLowRSquaredValues, 
                        noLightGaussLowRSquaredValues, 
                        midLightGaussLowRSquaredValues,
                        highLightGaussLowRSquaredValues):

    """
    Load dataset
    """
    [response_noLight, response_midLight, response_highLight, 
     sem_noLight, sem_midLight, sem_highLight] = load_ResponseData(responses_all,sem_all)
    
    numSoundResponsiveCells = responses_all.shape[1]
        
    miNoLight_allCells = np.zeros((numSoundResponsiveCells,))
    miMidLight_allCells = np.zeros((numSoundResponsiveCells,))
    miHighLight_allCells = np.zeros((numSoundResponsiveCells,))
                                    
    fitParamsSigmoidNoLight_allCells = np.zeros((4,numSoundResponsiveCells))
    fitParamsSigmoidMidLight_allCells = np.zeros((4,numSoundResponsiveCells))
    fitParamsSigmoidHighLight_allCells = np.zeros((4,numSoundResponsiveCells))
    rSquaredSigmoidNoLight_allCells = np.zeros((numSoundResponsiveCells,))
    rSquaredSigmoidMidLight_allCells = np.zeros((numSoundResponsiveCells,))
    rSquaredSigmoidHighLight_allCells = np.zeros((numSoundResponsiveCells,))

    fitParamsGaussNoLight_allCells = np.zeros((4,numSoundResponsiveCells))
    fitParamsGaussMidLight_allCells = np.zeros((4,numSoundResponsiveCells))
    fitParamsGaussHighLight_allCells = np.zeros((4,numSoundResponsiveCells))
    rSquaredGaussNoLight_allCells = np.zeros((numSoundResponsiveCells,))
    rSquaredGaussMidLight_allCells = np.zeros((numSoundResponsiveCells,))
    rSquaredGaussHighLight_allCells = np.zeros((numSoundResponsiveCells,))
    
    """
    Bounds for minimization
    """
    sigBnds = [(-np.inf,np.inf),(-np.inf,np.inf),(-np.inf,np.inf),(5,np.inf)]
    gaussBnds = [(-1.2,1.2),(0,90),(3.3,np.inf),(-np.inf,np.inf)]

    for i in range(numSoundResponsiveCells):
        
        """
        MI across conditions
        """
        miNoLight_allCells[i] = monotonicity_rawData(response_noLight[:,i])
        miMidLight_allCells[i] = monotonicity_rawData(response_midLight[:,i])
        miHighLight_allCells[i] = monotonicity_rawData(response_highLight[:,i])
        
        """
        Fitting procedure using scipy.optimize.fmin. 

        To do: take different starting points and take union
        """
        
        fitParamsSigmoidNoLight = scipy.optimize.minimize(sigmoid_fn, [-0.1,-0.4,40,7.5],
                                                            args=(response_noLight[:,i],sem_noLight[:,i]),
                                                            method = 'Powell', bounds = sigBnds,
                                                        options={'maxiter':30000, 'maxfev':30000,'ftol':0.01, 'disp':0}).x

        fitParamsSigmoidNoLight_allCells[:,i] = fitParamsSigmoidNoLight
        rSquaredSigmoidNoLight_allCells[i] = sigmoid_fn(fitParamsSigmoidNoLight_allCells[:,i],response_noLight[:,i],
                                                        sem_noLight[:,i])

        fitParamsGaussNoLight = scipy.optimize.minimize(gaussian_fn, [-0.57,35,5.8,0],
                                                args=(response_noLight[:,i],sem_noLight[:,i]),
                                                    method = 'Powell', bounds = gaussBnds,
                                                    options={'maxiter':30000, 'maxfev':30000,'ftol':0.01, 'disp':0}).x

        fitParamsGaussNoLight_allCells[:,i] = fitParamsGaussNoLight
        rSquaredGaussNoLight_allCells[i] = gaussian_fn(fitParamsGaussNoLight_allCells[:,i],response_noLight[:,i],
                                                       sem_noLight[:,i])

        fitParamsSigmoidMidLight = scipy.optimize.minimize(sigmoid_fn, [-0.1,-0.4,40,7.5],
                                                            args=(response_midLight[:,i],sem_midLight[:,i]),
                                                            method = 'Powell', bounds = sigBnds,
                                                        options={'maxiter':30000, 'maxfev':30000,'ftol':0.01, 'disp':0}).x
        
        fitParamsSigmoidMidLight_allCells[:,i] = fitParamsSigmoidMidLight
        rSquaredSigmoidMidLight_allCells[i] = sigmoid_fn(fitParamsSigmoidMidLight_allCells[:,i],response_midLight[:,i],
                                                         sem_midLight[:,i])

        fitParamsGaussMidLight = scipy.optimize.minimize(gaussian_fn, [-0.57,35,5.8,0],
                                                args=(response_midLight[:,i],sem_midLight[:,i]),
                                                    method = 'Powell', bounds = gaussBnds,
                                                    options={'maxiter':30000, 'maxfev':30000,'ftol':0.01, 'disp':0}).x

        fitParamsGaussMidLight_allCells[:,i] = fitParamsGaussMidLight
        rSquaredGaussMidLight_allCells[i] = gaussian_fn(fitParamsGaussMidLight_allCells[:,i],response_midLight[:,i],
                                                        sem_midLight[:,i])

        fitParamsSigmoidHighLight = scipy.optimize.minimize(sigmoid_fn, [-0.1,-0.4,40,7.5],
                                                            args=(response_highLight[:,i],sem_highLight[:,i]),
                                                            method = 'Powell', bounds = sigBnds,
                                                        options={'maxiter':30000, 'maxfev':30000,'ftol':0.01, 'disp':0}).x
        
        fitParamsSigmoidHighLight_allCells[:,i] = fitParamsSigmoidHighLight
        rSquaredSigmoidHighLight_allCells[i] = sigmoid_fn(fitParamsSigmoidHighLight_allCells[:,i],response_highLight[:,i],
                                                          sem_highLight[:,i])

        fitParamsGaussHighLight = scipy.optimize.minimize(gaussian_fn, [-0.57,35,5.8,0],
                                                args=(response_highLight[:,i],sem_highLight[:,i]),
                                                    method = 'Powell', bounds = gaussBnds,
                                                    options={'maxiter':30000, 'maxfev':30000,'ftol':0.01, 'disp':0}).x

        fitParamsGaussHighLight_allCells[:,i] = fitParamsGaussHighLight
        rSquaredGaussHighLight_allCells[i] = gaussian_fn(fitParamsGaussHighLight_allCells[:,i],response_highLight[:,i],
                                                         sem_highLight[:,i])
    
    """
    Nonmonotonic if MI  < 0.5, Monotonic if MI >= 0.8, maybe either if MI >=0.5 and MI < 0.8
    """
    miNoLightGauss = np.where(miNoLight_allCells < 0.3)[0]
    miNoLightSigmoid = np.where(miNoLight_allCells >= 0.7)[0]
    miNoLightIntersection = np.intersect1d(np.where(miNoLight_allCells < 0.7)[0],
                                           np.where(miNoLight_allCells >= 0.3)[0])
    miMidLightGauss = np.where(miMidLight_allCells < 0.3)[0]
    miMidLightSigmoid = np.where(miMidLight_allCells >= 0.7)[0] 
    miMidLightIntersection = np.intersect1d(np.where(miMidLight_allCells < 0.7)[0],
                                            np.where(miMidLight_allCells >= 0.3)[0])
    miHighLightGauss = np.where(miHighLight_allCells < 0.3)[0]
    miHighLightSigmoid = np.where(miHighLight_allCells >= 0.7)[0]
    miHighLightIntersection = np.intersect1d(np.where(miHighLight_allCells < 0.7)[0],
                                             np.where(miHighLight_allCells >= 0.3)[0])
    
    """
    Testing MI
    """
    #print("Cells to be fit using Gaussian function")
    #for j in miNoLightGauss:
    #    plt.plot(amplitude, response_noLight[:,j])
    #    plt.show()
    #    print(miNoLight_allCells[j])
    
    #print("Cells to be fit using Sigmoid function")
    #for j in miNoLightSigmoid:
    #    plt.plot(amplitude, response_noLight[:,j])
    #    plt.show()
    #    print(miNoLight_allCells[j])

    #Counting number of cells in each condition that are well fit by the sigmoid and the gaussian curves.    
    
    """
    No Light condition
    """
    idxNoLightSigmoidLowRSquared = np.intersect1d(np.where(rSquaredSigmoidNoLight_allCells<0.2)[0], 
                                                  miNoLightSigmoid)
    idxNoLightGaussLowRSquared = np.intersect1d(np.where(rSquaredGaussNoLight_allCells<0.2)[0],
                                                miNoLightGauss)
    
    idxNoLightIntersectionSigmoid = np.intersect1d(np.where(rSquaredSigmoidNoLight_allCells<0.2)[0],
                                                   miNoLightIntersection)
    idxNoLightIntersectionGauss = np.intersect1d(np.where(rSquaredGaussNoLight_allCells<0.2)[0],
                                                 miNoLightIntersection)
    
    noLightIntersect = np.intersect1d(idxNoLightIntersectionSigmoid,
                                      idxNoLightIntersectionGauss)
    IntersectFitBySigmoid = rSquaredSigmoidNoLight_allCells[noLightIntersect] < rSquaredGaussNoLight_allCells[noLightIntersect]
    noLightSigmoidUnique = idxNoLightIntersectionSigmoid[~np.isin(idxNoLightIntersectionSigmoid,
                                                                 noLightIntersect[~IntersectFitBySigmoid])]
    noLightGaussUnique = idxNoLightIntersectionGauss[~np.isin(idxNoLightIntersectionGauss,
                                                             noLightIntersect[IntersectFitBySigmoid])]
    
    noLightSigmoidComplete = np.concatenate((idxNoLightSigmoidLowRSquared,noLightSigmoidUnique))
    noLightGaussComplete = np.concatenate((idxNoLightGaussLowRSquared,noLightGaussUnique))
    #SigmoidToKeepCaseI = np.intersect1d(np.where(fitParamsSigmoidNoLight_allCells[3,:]>0)[0],
    #                               np.where(fitParamsSigmoidNoLight_allCells[1,:]>fitParamsSigmoidNoLight_allCells[0,:])[0])
    #SigmoidToKeepCaseII = np.intersect1d(np.where(fitParamsSigmoidNoLight_allCells[3,:]<0)[0],
    #                               np.where(fitParamsSigmoidNoLight_allCells[1,:]<fitParamsSigmoidNoLight_allCells[0,:])[0])
    #SigmoidToKeep = np.concatenate((SigmoidToKeepCaseI, SigmoidToKeepCaseII))
    #noLightSigmoidComplete_PosSlope = np.intersect1d(noLightSigmoidComplete, SigmoidToKeep)
    #noLightSigmoidLowRSquared.append(noLightSigmoidComplete_PosSlope)
    #noLightSigmoidLowRSquaredFitParams.append(fitParamsSigmoidNoLight_allCells[:,noLightSigmoidComplete_PosSlope])
    noLightSigmoidLowRSquared.append(noLightSigmoidComplete)
    noLightSigmoidLowRSquaredFitParams.append(fitParamsSigmoidNoLight_allCells[:,noLightSigmoidComplete])
    noLightSigmoidLowRSquaredValues.append(rSquaredSigmoidNoLight_allCells[noLightSigmoidComplete])
    
    GaussianToKeep = np.intersect1d(noLightGaussComplete[fitParamsGaussNoLight_allCells[1,noLightGaussComplete]<89],
                                    noLightGaussComplete[fitParamsGaussNoLight_allCells[1,noLightGaussComplete]>1])
    noLightGaussLowRSquared.append(GaussianToKeep)
    noLightGaussLowRSquaredFitParams.append(fitParamsGaussNoLight_allCells[:,GaussianToKeep])
    noLightGaussLowRSquaredValues.append(rSquaredGaussNoLight_allCells[GaussianToKeep])
    
    """
    for j in GaussianToKeep: 
        plt.figure()
        [amp0,mean0,std0,shift] = fitParamsGaussNoLight_allCells[:,j]
        y = shift + amp0*np.exp(-(amplitude-mean0)**2/(2*std0**2))
        plt.plot(amplitude, response_noLight[:,j],'k')
        plt.plot(amplitude, y,'k--') 
        plt.show()
        print(j, miNoLight_allCells[j], fitParamsGaussNoLight_allCells[:,j])
        
    for j in noLightSigmoidComplete:
        plt.figure()
        [y0, ymax, x0, delta_x] = fitParamsSigmoidNoLight_allCells[:,j]
        y = y0 + (ymax-y0)/(1+np.exp((x0-amplitude)/delta_x)); 
        plt.plot(amplitude, response_noLight[:,j],'k')
        plt.plot(amplitude, y,'k--') 
        plt.show()
        print(j, miNoLight_allCells[j], fitParamsSigmoidNoLight_allCells[:,j])    
    """
    
    """
    Mid Light condition
    """
    idxMidLightSigmoidLowRSquared = np.intersect1d(np.where(rSquaredSigmoidMidLight_allCells<0.2)[0], 
                                                  miMidLightSigmoid)
    idxMidLightGaussLowRSquared = np.intersect1d(np.where(rSquaredGaussMidLight_allCells<0.2)[0],
                                                miMidLightGauss)
    
    idxMidLightIntersectionSigmoid = np.intersect1d(np.where(rSquaredSigmoidMidLight_allCells<0.2)[0],
                                                   miMidLightIntersection)
    idxMidLightIntersectionGauss = np.intersect1d(np.where(rSquaredGaussMidLight_allCells<0.2)[0],
                                                 miMidLightIntersection)
    
    midLightIntersect = np.intersect1d(idxMidLightIntersectionSigmoid,
                                      idxMidLightIntersectionGauss)
    IntersectFitBySigmoid = rSquaredSigmoidMidLight_allCells[midLightIntersect] < rSquaredGaussMidLight_allCells[midLightIntersect]
    midLightSigmoidUnique = idxMidLightIntersectionSigmoid[~np.isin(idxMidLightIntersectionSigmoid,
                                                                 midLightIntersect[~IntersectFitBySigmoid])]
    midLightGaussUnique = idxMidLightIntersectionGauss[~np.isin(idxMidLightIntersectionGauss,
                                                             midLightIntersect[IntersectFitBySigmoid])]
    
    midLightSigmoidComplete = np.concatenate((idxMidLightSigmoidLowRSquared,midLightSigmoidUnique))
    midLightGaussComplete = np.concatenate((idxMidLightGaussLowRSquared,midLightGaussUnique))
    
    midLightSigmoidLowRSquared.append(midLightSigmoidComplete)
    midLightSigmoidLowRSquaredFitParams.append(fitParamsSigmoidMidLight_allCells[:,midLightSigmoidComplete])
    midLightSigmoidLowRSquaredValues.append(rSquaredSigmoidMidLight_allCells[midLightSigmoidComplete])                                    
    GaussianToKeep = np.intersect1d(midLightGaussComplete[fitParamsGaussMidLight_allCells[1,midLightGaussComplete]<89],
                                    midLightGaussComplete[fitParamsGaussMidLight_allCells[1,midLightGaussComplete]>1])
    midLightGaussLowRSquared.append(GaussianToKeep)
    midLightGaussLowRSquaredFitParams.append(fitParamsGaussMidLight_allCells[:,GaussianToKeep])
    midLightGaussLowRSquaredValues.append(rSquaredGaussMidLight_allCells[GaussianToKeep])
    
    """
    for j in midLightGaussComplete:        
        [amp0,mean0,std0,shift] = fitParamsGaussMidLight_allCells[:,j]
        y = shift + amp0*np.exp(-(amplitude-mean0)**2/(2*std0**2))
        plt.plot(amplitude, response_midLight[:,j],'k')
        plt.plot(amplitude, y,'k--') 
        plt.show()
        print(miMidLight_allCells[j])
        
    for j in midLightSigmoidComplete:
        [y0, ymax, x0, delta_x] = fitParamsSigmoidMidLight_allCells[:,j]
        y = y0 + (ymax-y0)/(1+np.exp((x0-amplitude)/delta_x)); 
        plt.plot(amplitude, response_midLight[:,j],'k')
        plt.plot(amplitude, y,'k--') 
        plt.show()
        print(miMidLight_allCells[j]) 
    """
    
    """
    High Light condition
    """
    idxHighLightSigmoidLowRSquared = np.intersect1d(np.where(rSquaredSigmoidHighLight_allCells<0.2)[0], 
                                                  miHighLightSigmoid)
    idxHighLightGaussLowRSquared = np.intersect1d(np.where(rSquaredGaussHighLight_allCells<0.2)[0],
                                                miHighLightGauss)
    
    idxHighLightIntersectionSigmoid = np.intersect1d(np.where(rSquaredSigmoidHighLight_allCells<0.2)[0],
                                                   miHighLightIntersection)
    idxHighLightIntersectionGauss = np.intersect1d(np.where(rSquaredGaussHighLight_allCells<0.2)[0],
                                                 miHighLightIntersection)
    
    highLightIntersect = np.intersect1d(idxHighLightIntersectionSigmoid,
                                      idxHighLightIntersectionGauss)
    IntersectFitBySigmoid = rSquaredSigmoidHighLight_allCells[highLightIntersect] < rSquaredGaussHighLight_allCells[highLightIntersect]
    highLightSigmoidUnique = idxHighLightIntersectionSigmoid[~np.isin(idxHighLightIntersectionSigmoid,
                                                                 highLightIntersect[~IntersectFitBySigmoid])]
    highLightGaussUnique = idxHighLightIntersectionGauss[~np.isin(idxHighLightIntersectionGauss,
                                                             highLightIntersect[IntersectFitBySigmoid])]
    
    highLightSigmoidComplete = np.concatenate((idxHighLightSigmoidLowRSquared,highLightSigmoidUnique))
    highLightGaussComplete = np.concatenate((idxHighLightGaussLowRSquared,highLightGaussUnique))
    
    highLightSigmoidLowRSquared.append(highLightSigmoidComplete)
    highLightSigmoidLowRSquaredFitParams.append(fitParamsSigmoidHighLight_allCells[:,highLightSigmoidComplete])
    highLightSigmoidLowRSquaredValues.append(rSquaredSigmoidHighLight_allCells[highLightSigmoidComplete])                   
    
    GaussianToKeep = np.intersect1d(highLightGaussComplete[fitParamsGaussHighLight_allCells[1,highLightGaussComplete]<89],
                                    highLightGaussComplete[fitParamsGaussHighLight_allCells[1,highLightGaussComplete]>1])
    highLightGaussLowRSquared.append(GaussianToKeep)
    highLightGaussLowRSquaredFitParams.append(fitParamsGaussHighLight_allCells[:,GaussianToKeep])
    highLightGaussLowRSquaredValues.append(rSquaredGaussHighLight_allCells[GaussianToKeep])  
    
    """
    for j in highLightGaussComplete:        
        [amp0,mean0,std0,shift] = fitParamsGaussHighLight_allCells[:,j]
        y = shift + amp0*np.exp(-(amplitude-mean0)**2/(2*std0**2))
        plt.plot(amplitude, response_highLight[:,j],'k')
        plt.plot(amplitude, y,'k--') 
        plt.show()
        print(miHighLight_allCells[j])
        
    for j in highLightSigmoidComplete:
        [y0, ymax, x0, delta_x] = fitParamsSigmoidHighLight_allCells[:,j]
        y = y0 + (ymax-y0)/(1+np.exp((x0-amplitude)/delta_x)); 
        plt.plot(amplitude, response_highLight[:,j],'k')
        plt.plot(amplitude, y,'k--') 
        plt.show()
        print(miHighLight_allCells[j])
    """
    
    return [noLightSigmoidLowRSquared, midLightSigmoidLowRSquared, highLightSigmoidLowRSquared, 
            noLightGaussLowRSquared, midLightGaussLowRSquared, highLightGaussLowRSquared,
            noLightSigmoidLowRSquaredFitParams, midLightSigmoidLowRSquaredFitParams, 
            highLightSigmoidLowRSquaredFitParams, noLightGaussLowRSquaredFitParams, 
            midLightGaussLowRSquaredFitParams, highLightGaussLowRSquaredFitParams,
            noLightSigmoidLowRSquaredValues, midLightSigmoidLowRSquaredValues, 
            highLightSigmoidLowRSquaredValues, noLightGaussLowRSquaredValues, 
            midLightGaussLowRSquaredValues, highLightGaussLowRSquaredValues]
    