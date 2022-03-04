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
from utils import *

def fittingStatisticsGauss(responseMatrix,
                           num_cells,
                          lightCondition,
                          listGaussAllConds,
                          listGaussAllCondsStart,
                          dictGaussFitParams,
                          dictGaussIndices,
                          dictGaussRSquaredValues,
                          listSigmoidAllConds,
                          listSigmoidAllCondsStart,
                          dictSigmoidFitParams,
                          dictSigmoidIndices,
                          dictSigmoidRSquaredValues):
    
    gaussParametersArray = np.zeros((4,1))
    interpolatedErrorGauss = np.zeros((1))
    gaussDatasetIndex = np.zeros((1))
    gaussCellIndex = np.zeros((1))
    sigmoidParametersArray = np.zeros((4,1))
    sigmoidY0Array = np.zeros((1))
    sigmoidY90Array = np.zeros((1))
    interpolatedErrorSigmoid = np.zeros((1))
    sigmoidDatasetIndex = np.zeros((1))
    sigmoidCellIndex = np.zeros((1))
    cntGaussPositive = 0 # variable that counts how many cells fit using gauss have positive MI
    cntSigmoidPositive = 0 # variable that counts how many cells fit using sigmoid have positive MI
    
    if lightCondition == 'no':
        indices = [0,3,6,9,12,15,18]
    elif lightCondition == 'mid':
        indices = [1,4,7,10,13,16,19]
    else:
        indices = [2,5,8,11,14,17,20]
        
    x_sampled = np.arange(0,90,2)    
        
    for dataset in range(num_cells): 
        for posCells in range(sum(listSigmoidAllConds[dataset,:]>-1)):
            if listSigmoidAllConds[dataset,posCells] > -1:
                cell = int(listSigmoidAllConds[dataset,posCells])
                startPt = int(listSigmoidAllCondsStart[dataset,posCells])
                bestParameters = dictSigmoidFitParams[startPt,dataset][0][:,dictSigmoidIndices[startPt,dataset][0]==cell] 
                mcFaddenRSquared = dictSigmoidRSquaredValues[startPt,dataset][0][dictSigmoidIndices[startPt,dataset][0]==cell]
                avgResponse = responseMatrix['Ind_CellSelectionIncludingSEM'][dataset,0][indices,cell]
                semResponse = responseMatrix['Ind_CellSelectionIncludingSEM'][dataset,1][indices,cell]
                ifSigSoundPositive,newBestParameters = check_soundPositive(parameters=bestParameters,
                                                                           ifGauss=0)
                newBestParameters[:2] = newBestParameters[:2]*np.max(np.abs(avgResponse))
                sigmoidParametersArray = np.append(sigmoidParametersArray,newBestParameters,axis=1)
                cntSigmoidPositive += int(ifSigSoundPositive)
                interpolated_err, _, _ = interpolatedSigmoid_Rsquared(params=newBestParameters,response=avgResponse)
                interpolatedErrorSigmoid = np.append(interpolatedErrorSigmoid,interpolated_err)   
                sigmoidDatasetIndex = np.append(sigmoidDatasetIndex,dataset)
                sigmoidCellIndex = np.append(sigmoidCellIndex,cell)
                
                plt.figure()
                plt.errorbar(amplitude, avgResponse, yerr=semResponse)
                [y0, ymax, x0, delta_x] = newBestParameters
                y = y0 + (ymax-y0)/(1+np.exp((x0-x_sampled)/delta_x))
                plt.plot(x_sampled, y,'k--') 
                y = y0 + (ymax-y0)/(1+np.exp((x0-amplitude)/delta_x))
                plt.plot(amplitude, y,'b--') 
                plt.xlabel('Intensity',fontsize=15)
                plt.ylabel('Response',fontsize=15)
                plt.xticks(fontsize=15)
                plt.yticks(fontsize=15)
                plt.tight_layout()
                print('params',newBestParameters, 
                      'MI',monotonicity_rawData(avgResponse), 
                      'Sound Positive?', ifSigSoundPositive, 
                      'McFadden R-squared', mcFaddenRSquared,
                      'Interpolated error', interpolated_err)
                
                [y0, ymax, x0, delta_x] = newBestParameters
                sigmoidY0Array = np.append(sigmoidY0Array,y0 + (ymax-y0)/(1+np.exp((x0-0)/delta_x)))
                sigmoidY90Array = np.append(sigmoidY90Array,y0 + (ymax-y0)/(1+np.exp((x0-90)/delta_x)))
                
        for posCells in range(sum(listGaussAllConds[dataset,:]>-1)):
            if listGaussAllConds[dataset,posCells] > -1:
                cell = int(listGaussAllConds[dataset,posCells])
                startPt = int(listGaussAllCondsStart[dataset,posCells])
                bestParameters = dictGaussFitParams[startPt,dataset][0][:,dictGaussIndices[startPt,dataset][0]==cell]          
                mcFaddenRSquared = dictGaussRSquaredValues[startPt,dataset][0][dictGaussIndices[startPt,dataset][0]==cell]       
                avgResponse = responseMatrix['Ind_CellSelectionIncludingSEM'][dataset,0][indices,cell]
                semResponse = responseMatrix['Ind_CellSelectionIncludingSEM'][dataset,1][indices,cell]
                ifGaussSoundPositive, _ = check_soundPositive(parameters=bestParameters,
                                                              ifGauss=1)
                bestParameters[[0,3]] = bestParameters[[0,3]]*np.max(np.abs(avgResponse))
                if 10 <= bestParameters[1] <= 80:
                    gaussParametersArray = np.append(gaussParametersArray,bestParameters,axis=1)
                    cntGaussPositive += int(ifGaussSoundPositive)
                    interpolated_err, interpolated_amp, interpolated_response = interpolatedGauss_Rsquared(params=bestParameters, response=avgResponse)
                    interpolatedErrorGauss = np.append(interpolatedErrorGauss,interpolated_err)
                    gaussDatasetIndex = np.append(gaussDatasetIndex,dataset)
                    gaussCellIndex = np.append(gaussCellIndex,cell)
                    
                    if interpolated_err:
                        plt.errorbar(amplitude, avgResponse, yerr=semResponse)
                        plt.plot(interpolated_amp, interpolated_response,'o')
                        [amp0,mean0,std0,shift] = bestParameters  
                        y = shift+amp0*np.exp(-(x_sampled-mean0)**2/(2*std0**2))
                        plt.plot(x_sampled, y,'k--') 
                        y = shift+amp0*np.exp(-(amplitude-mean0)**2/(2*std0**2))
                        plt.plot(amplitude, y,'r--') 
                        plt.xlabel('Intensity',fontsize=15)
                        plt.ylabel('Response',fontsize=15)
                        plt.xticks(fontsize=15)
                        plt.yticks(fontsize=15)
                        plt.tight_layout()
                        print('params',bestParameters, 
                              'MI',monotonicity_rawData(avgResponse), 
                              'Sound Positive?', ifGaussSoundPositive, 
                              'McFadden R-squared', mcFaddenRSquared,
                              'Interpolated error', interpolated_err)
                        #pdb.set_trace()
                        plt.show()  
                        
                            
                else:
                    plt.figure()
                    plt.errorbar(amplitude, avgResponse, yerr=semResponse)
                    [amp0,mean0,std0,shift] = bestParameters  
                    y = shift+amp0*np.exp(-(x_sampled-mean0)**2/(2*std0**2))
                    plt.plot(x_sampled, y,'k--') 
                    y = shift+amp0*np.exp(-(amplitude-mean0)**2/(2*std0**2))
                    plt.plot(amplitude, y,'r--') 
                    sigmoidFitParameters = moveToSigmoid(avgResponse,semResponse)
                    print(sigmoidFitParameters)
                    ifSigSoundPositive,newSigmoidParameters = check_soundPositive(parameters=sigmoidFitParameters,
                                                                                   ifGauss=0)
                    newSigmoidParameters[:2] = newSigmoidParameters[:2]*np.max(np.abs(avgResponse))
                    cntSigmoidPositive += int(ifSigSoundPositive)
                    sigmoidParametersArray = np.append(sigmoidParametersArray,newSigmoidParameters[:,np.newaxis],axis=1)
                    interpolated_err, _, _ = interpolatedSigmoid_Rsquared(params=newSigmoidParameters, 
                                                                          response=avgResponse)
                    interpolatedErrorSigmoid = np.append(interpolatedErrorSigmoid,interpolated_err)
                    sigmoidDatasetIndex = np.append(sigmoidDatasetIndex,dataset)
                    sigmoidCellIndex = np.append(sigmoidCellIndex,cell)
                    [y0, ymax, x0, delta_x] = newSigmoidParameters                    
                    sigmoidY0Array = np.append(sigmoidY0Array,y0 + (ymax-y0)/(1+np.exp((x0-0)/delta_x)))
                    sigmoidY90Array = np.append(sigmoidY90Array,y0 + (ymax-y0)/(1+np.exp((x0-90)/delta_x)))
                    y = y0 + (ymax-y0)/(1+np.exp((x0-x_sampled)/delta_x))
                    plt.plot(x_sampled, y,'k--') 
                    y = y0 + (ymax-y0)/(1+np.exp((x0-amplitude)/delta_x))
                    plt.plot(amplitude, y,'b--') 
                    #print('params',newSigmoidParameters, 
                    #  'MI',monotonicity_rawData(avgResponse), 
                    #  'Sound Positive?', ifSigSoundPositive, 
                    #  'McFadden R-squared', mcFaddenRSquared,
                    #  'Interpolated error', interpolated_err,
                    #  'Value at 0 dB', y0 + (ymax-y0)/(1+np.exp((x0-0)/delta_x)),
                    #  'Value at 90 dB', y0 + (ymax-y0)/(1+np.exp((x0-90)/delta_x)))
                    plt.show()
                    
                
    return gaussParametersArray[:,1:], cntGaussPositive, interpolatedErrorGauss[1:], gaussDatasetIndex[1:], gaussCellIndex[1:], sigmoidParametersArray[:,1:], sigmoidY0Array[1:], sigmoidY90Array[1:], cntSigmoidPositive, interpolatedErrorSigmoid[1:], sigmoidDatasetIndex[1:], sigmoidCellIndex[1:]


def moveToSigmoid(response_Light,
                  sem_Light):
    
    startingValues = [[0.4,-0.4,11,-5],[0.4,-0.4,11,10],[0.4,-1,11,-5],[0.4,-1,11,10],[-0.4,-0.6,40,-5],[-0.4,-0.6,40,10],
                      [-0.8,0.2,70,-5],[-0.8,0.2,70,10],[0.2,0.8,70,-5],[0.2,0.8,70,7.5],[-0.6,0.6,30,-5],[-0.6,0.6,30,7.5],
                      [-0.1,-0.5,40,-5],[-0.1,-0.5,40,7.5],[0,0.5,70,-5],[0,0.5,70,7.5]]    
    
    fitParamsSigmoid = np.empty((4,len(startingValues)))
    rSquaredSigmoid = np.empty(len(startingValues))
    sigBndsPos = [(-np.inf,np.inf),(-np.inf,np.inf),(-np.inf,np.inf),(5,np.inf)]
    sigBndsNeg = [(-np.inf,np.inf),(-np.inf,np.inf),(-np.inf,np.inf),(-np.inf,-5)]
    for iStart in range(len(startingValues)):
        if startingValues[iStart][3]>0:
            fitParamsSigmoid[:,iStart] = scipy.optimize.minimize(sigmoid_fn, startingValues[iStart],
                                                           args=(response_Light/np.max(np.abs(response_Light)),
                                                                 sem_Light/np.max(np.abs(response_Light))),
                                                           method = 'Powell', bounds = sigBndsPos,
                                                           options={'maxiter':30000, 'maxfev':30000,'ftol':0.01, 'disp':0}).x

            rSquaredSigmoid[iStart] = sigmoid_fn(fitParamsSigmoid[:,iStart],response_Light/np.max(np.abs(response_Light)), sem_Light/np.max(np.abs(response_Light)))
            
        if startingValues[iStart][3]<0:
            fitParamsSigmoid[:,iStart] = scipy.optimize.minimize(sigmoid_fn, startingValues[iStart],
                                                           args=(response_Light/np.max(np.abs(response_Light)),
                                                                 sem_Light/np.max(np.abs(response_Light))),
                                                           method = 'Powell', bounds = sigBndsNeg,
                                                           options={'maxiter':30000, 'maxfev':30000,'ftol':0.01, 'disp':0}).x

            rSquaredSigmoid[iStart] = sigmoid_fn(fitParamsSigmoid[:,iStart],response_Light/np.max(np.abs(response_Light)), sem_Light/np.max(np.abs(response_Light)))

    return fitParamsSigmoid[:,np.argmin(rSquaredSigmoid)]
    
    
