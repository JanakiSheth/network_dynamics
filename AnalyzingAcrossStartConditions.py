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

"""
Combining different starting point conditions for all the datasets across different light conditions

Comparing R-squared across starting conditions and across sigmoid and gaussian cases for the overlapping cells. Then plotting the best fit against an errorbar plot which uses average and sem. There are some bad fits which have passed the mi test and r-squared test. The reason they are bad is because of high sem. One way of overcoming this is to eyeball how the fit looks. 
"""

def combiningStartPts(responseFunctionMat,
                    num_cells, 
                    dictNoLightSigmoidLowRSquared, 
                    dictMidLightSigmoidLowRSquared, 
                    dictHighLightSigmoidLowRSquared, 
                    dictNoLightGaussLowRSquared, 
                    dictMidLightGaussLowRSquared, 
                    dictHighLightGaussLowRSquared,
                    dictNoLightSigmoidLowRSquaredFitParams, 
                    dictMidLightSigmoidLowRSquaredFitParams, 
                    dictHighLightSigmoidLowRSquaredFitParams, 
                    dictNoLightGaussLowRSquaredFitParams,
                    dictMidLightGaussLowRSquaredFitParams, 
                    dictHighLightGaussLowRSquaredFitParams,
                    dictNoLightSigmoidLowRSquaredValues, 
                    dictMidLightSigmoidLowRSquaredValues, 
                    dictHighLightSigmoidLowRSquaredValues, 
                    dictNoLightGaussLowRSquaredValues, 
                    dictMidLightGaussLowRSquaredValues,
                    dictHighLightGaussLowRSquaredValues,
                    num_startingPts = 7):

    listNoLightGaussAllConds = np.ones((num_cells,80))*-1
    listMidLightGaussAllConds = np.ones((num_cells,80))*-1
    listHighLightGaussAllConds = np.ones((num_cells,80))*-1
    listNoLightGaussAllCondsStart = np.ones((num_cells,80))*-1
    listMidLightGaussAllCondsStart = np.ones((num_cells,80))*-1
    listHighLightGaussAllCondsStart = np.ones((num_cells,80))*-1

    listNoLightSigmoidAllConds = np.ones((num_cells,80))*-1
    listMidLightSigmoidAllConds = np.ones((num_cells,80))*-1
    listHighLightSigmoidAllConds = np.ones((num_cells,80))*-1
    listNoLightSigmoidAllCondsStart = np.ones((num_cells,80))*-1
    listMidLightSigmoidAllCondsStart = np.ones((num_cells,80))*-1
    listHighLightSigmoidAllCondsStart = np.ones((num_cells,80))*-1

    for ii in range(num_cells):
        listNoLightGaussAllCells = np.array([])
        listMidLightGaussAllCells = np.array([])
        listHighLightGaussAllCells = np.array([])
        listNoLightGaussAllCellsRValue = np.array([])
        listMidLightGaussAllCellsRValue = np.array([])
        listHighLightGaussAllCellsRValue = np.array([])
        listNoLightGaussStartingPt = np.array([])
        listMidLightGaussStartingPt = np.array([])
        listHighLightGaussStartingPt = np.array([])

        """
        Comparing R-squared value across different starting conditions and storing starting point corresponding
        to the smallest R-squared value
        """

        for jj in range(num_startingPts):
            listNoLightGaussAllCells = np.append(listNoLightGaussAllCells,
                                                 dictNoLightGaussLowRSquared[jj,ii])
            listNoLightGaussAllCellsRValue = np.append(listNoLightGaussAllCellsRValue,
                                                 dictNoLightGaussLowRSquaredValues[jj,ii])
            listNoLightGaussStartingPt = np.append(listNoLightGaussStartingPt,
                                              np.ones((np.shape(dictNoLightGaussLowRSquared[jj,ii])[1],))*jj)
            listMidLightGaussAllCells = np.append(listMidLightGaussAllCells,
                                                  dictMidLightGaussLowRSquared[jj,ii])
            listMidLightGaussAllCellsRValue = np.append(listMidLightGaussAllCellsRValue,
                                                 dictMidLightGaussLowRSquaredValues[jj,ii])
            listMidLightGaussStartingPt = np.append(listMidLightGaussStartingPt,
                                              np.ones((np.shape(dictMidLightGaussLowRSquared[jj,ii])[1],))*jj)
            listHighLightGaussAllCells = np.append(listHighLightGaussAllCells,
                                                   dictHighLightGaussLowRSquared[jj,ii])
            listHighLightGaussAllCellsRValue = np.append(listHighLightGaussAllCellsRValue,
                                                 dictHighLightGaussLowRSquaredValues[jj,ii])
            listHighLightGaussStartingPt = np.append(listHighLightGaussStartingPt,
                                              np.ones((np.shape(dictHighLightGaussLowRSquared[jj,ii])[1],))*jj)

        if len(listNoLightGaussAllCells):
            listNoLightGaussAllConds[ii,:len(np.unique(listNoLightGaussAllCells))] = np.unique(listNoLightGaussAllCells) 
         
        kk = 0
        while listNoLightGaussAllConds[ii,kk]>-1:
            listNoLightGaussAllCondsStart[ii,kk] = listNoLightGaussStartingPt[np.multiply(
                listNoLightGaussAllCellsRValue==np.min(listNoLightGaussAllCellsRValue [listNoLightGaussAllCells==
                                                                                       listNoLightGaussAllConds[ii,kk]]),
                listNoLightGaussAllCells==listNoLightGaussAllConds[ii,kk])][0]
            kk+=1
 

        if len(listMidLightGaussAllCells):
            listMidLightGaussAllConds[ii,:len(np.unique(listMidLightGaussAllCells))] = np.unique(listMidLightGaussAllCells)
            
        kk = 0
        while listMidLightGaussAllConds[ii,kk]>-1:
            listMidLightGaussAllCondsStart[ii,kk] = listMidLightGaussStartingPt[np.multiply(
                listMidLightGaussAllCellsRValue==np.min(listMidLightGaussAllCellsRValue[listMidLightGaussAllCells==
                                                                                        listMidLightGaussAllConds[ii,kk]]),
                listMidLightGaussAllCells==listMidLightGaussAllConds[ii,kk])][0]
            kk+=1

        if len(listHighLightGaussAllCells):
            listHighLightGaussAllConds[ii,:len(np.unique(listHighLightGaussAllCells))] = np.unique(listHighLightGaussAllCells)
            
        kk = 0
        while listHighLightGaussAllConds[ii,kk]>-1:
            listHighLightGaussAllCondsStart[ii,kk] = listHighLightGaussStartingPt[np.multiply(
                listHighLightGaussAllCellsRValue==np.min(listHighLightGaussAllCellsRValue[listHighLightGaussAllCells==
                                                                                          listHighLightGaussAllConds[ii,kk]]),
                listHighLightGaussAllCells==listHighLightGaussAllConds[ii,kk])][0]
            kk+=1

        listNoLightSigmoidAllCells = np.array([])
        listMidLightSigmoidAllCells = np.array([])
        listHighLightSigmoidAllCells = np.array([])
        listNoLightSigmoidAllCellsRValue = np.array([])
        listMidLightSigmoidAllCellsRValue = np.array([])
        listHighLightSigmoidAllCellsRValue = np.array([])
        listNoLightSigmoidStartingPt = np.array([])
        listMidLightSigmoidStartingPt = np.array([])
        listHighLightSigmoidStartingPt = np.array([])

        for jj in range(num_startingPts):
            listNoLightSigmoidAllCells = np.append(listNoLightSigmoidAllCells,
                                                   dictNoLightSigmoidLowRSquared[jj,ii])        
            listNoLightSigmoidAllCellsRValue = np.append(listNoLightSigmoidAllCellsRValue,
                                                 dictNoLightSigmoidLowRSquaredValues[jj,ii])
            listNoLightSigmoidStartingPt = np.append(listNoLightSigmoidStartingPt,
                                              np.ones((np.shape(dictNoLightSigmoidLowRSquared[jj,ii])[1],))*jj)
            
            listMidLightSigmoidAllCells = np.append(listMidLightSigmoidAllCells,
                                                    dictMidLightSigmoidLowRSquared[jj,ii])
            listMidLightSigmoidAllCellsRValue = np.append(listMidLightSigmoidAllCellsRValue,
                                                 dictMidLightSigmoidLowRSquaredValues[jj,ii])
            listMidLightSigmoidStartingPt = np.append(listMidLightSigmoidStartingPt,
                                              np.ones((np.shape(dictMidLightSigmoidLowRSquared[jj,ii])[1],))*jj)
            
            listHighLightSigmoidAllCells = np.append(listHighLightSigmoidAllCells,
                                                     dictHighLightSigmoidLowRSquared[jj,ii])
            listHighLightSigmoidAllCellsRValue = np.append(listHighLightSigmoidAllCellsRValue,
                                                 dictHighLightSigmoidLowRSquaredValues[jj,ii])
            listHighLightSigmoidStartingPt = np.append(listHighLightSigmoidStartingPt,
                                              np.ones((np.shape(dictHighLightSigmoidLowRSquared[jj,ii])[1],))*jj)
        
        if len(listNoLightSigmoidAllCells):
            listNoLightSigmoidAllConds[ii,:len(np.unique(listNoLightSigmoidAllCells))] = np.unique(listNoLightSigmoidAllCells)      
        kk = 0
        while listNoLightSigmoidAllConds[ii,kk]>-1:
            listNoLightSigmoidAllCondsStart[ii,kk] = listNoLightSigmoidStartingPt[np.multiply(
                listNoLightSigmoidAllCellsRValue==np.min(
                    listNoLightSigmoidAllCellsRValue[listNoLightSigmoidAllCells==listNoLightSigmoidAllConds[ii,kk]]),
                listNoLightSigmoidAllCells==listNoLightSigmoidAllConds[ii,kk])][0]
            kk+=1

        if len(listMidLightSigmoidAllCells):
            listMidLightSigmoidAllConds[ii,:len(np.unique(listMidLightSigmoidAllCells))] = np.unique(listMidLightSigmoidAllCells)  
            
        kk = 0
        while listMidLightSigmoidAllConds[ii,kk]>-1:
            listMidLightSigmoidAllCondsStart[ii,kk] = listMidLightSigmoidStartingPt[np.multiply(
                                                                                     listMidLightSigmoidAllCellsRValue==np.min(
                                                                                         listMidLightSigmoidAllCellsRValue[listMidLightSigmoidAllCells==listMidLightSigmoidAllConds[ii,kk]]),
                listMidLightSigmoidAllCells==listMidLightSigmoidAllConds[ii,kk])][0]
            kk+=1

        if len(listHighLightSigmoidAllCells):
            listHighLightSigmoidAllConds[ii,:len(np.unique(listHighLightSigmoidAllCells))] = np.unique(listHighLightSigmoidAllCells)

        kk = 0
        while listHighLightSigmoidAllConds[ii,kk]>-1:
            listHighLightSigmoidAllCondsStart[ii,kk] = listHighLightSigmoidStartingPt[np.multiply(
                listHighLightSigmoidAllCellsRValue==np.min(
                    listHighLightSigmoidAllCellsRValue[listHighLightSigmoidAllCells==listHighLightSigmoidAllConds[ii,kk]]),
                listHighLightSigmoidAllCells==listHighLightSigmoidAllConds[ii,kk])][0]
            kk+=1

        """
        Finding common cells between sigmoid and gaussian and removing the low R-squared value cell from
        respective type
        """    
        cellsCommon = np.intersect1d(listNoLightSigmoidAllConds[ii,:][listNoLightSigmoidAllConds[ii,:]>-1],
                                   listNoLightGaussAllConds[ii,:][listNoLightGaussAllConds[ii,:]>-1])

        if len(cellsCommon):
            for idxCommon in cellsCommon:

                GaussIdx = dictNoLightGaussLowRSquared[int(listNoLightGaussAllCondsStart[ii,:]
                                                        [listNoLightGaussAllConds[ii,:]==idxCommon]),ii][0]
                GaussRValue = dictNoLightGaussLowRSquaredValues[int(listNoLightGaussAllCondsStart[ii,:]
                                                        [listNoLightGaussAllConds[ii,:]==idxCommon]),ii][0][GaussIdx==idxCommon]

                SigmoidIdx = dictNoLightSigmoidLowRSquared[int(listNoLightSigmoidAllCondsStart[ii,:]
                                                        [listNoLightSigmoidAllConds[ii,:]==idxCommon]),ii][0]
                SigmoidRValue = dictNoLightSigmoidLowRSquaredValues[int(listNoLightSigmoidAllCondsStart[ii,:]
                                                        [listNoLightSigmoidAllConds[ii,:]==idxCommon]),ii][0][SigmoidIdx==idxCommon]
                if GaussRValue < SigmoidRValue:
                    listNoLightGaussAllConds[ii,:][listNoLightGaussAllConds[ii,:]==idxCommon] = -1
                    listNoLightGaussAllCondsStart[ii,:][listNoLightGaussAllConds[ii,:]==idxCommon] = -1
                else:
                    listNoLightSigmoidAllConds[ii,:][listNoLightSigmoidAllConds[ii,:]==idxCommon] = -1
                    listNoLightSigmoidAllCondsStart[ii,:][listNoLightSigmoidAllConds[ii,:]==idxCommon] = -1                  
            sortIdxs = np.argsort(listNoLightGaussAllConds[ii,:])[::-1]    
            listNoLightGaussAllConds[ii,:] = listNoLightGaussAllConds[ii,sortIdxs]
            listNoLightGaussAllCondsStart[ii,:] = listNoLightGaussAllCondsStart[ii,sortIdxs]
            sortIdxs = np.argsort(listNoLightSigmoidAllConds[ii,:])[::-1]            
            listNoLightSigmoidAllConds[ii,:] = listNoLightSigmoidAllConds[ii,sortIdxs]
            listNoLightSigmoidAllCondsStart[ii,:] = listNoLightSigmoidAllCondsStart[ii,sortIdxs]


        cellsCommon = np.intersect1d(listMidLightSigmoidAllConds[ii,:][listMidLightSigmoidAllConds[ii,:]>-1],
                                   listMidLightGaussAllConds[ii,:][listMidLightGaussAllConds[ii,:]>-1])
        if len(cellsCommon):
            for idxCommon in cellsCommon:

                GaussIdx = dictMidLightGaussLowRSquared[int(listMidLightGaussAllCondsStart[ii,:]
                                                        [listMidLightGaussAllConds[ii,:]==idxCommon]),ii][0]
                GaussRValue = dictMidLightGaussLowRSquaredValues[int(listMidLightGaussAllCondsStart[ii,:]
                                                        [listMidLightGaussAllConds[ii,:]==idxCommon]),ii][0][GaussIdx==idxCommon]

                SigmoidIdx = dictMidLightSigmoidLowRSquared[int(listMidLightSigmoidAllCondsStart[ii,:]
                                                        [listMidLightSigmoidAllConds[ii,:]==idxCommon]),ii][0]
                SigmoidRValue = dictMidLightSigmoidLowRSquaredValues[int(listMidLightSigmoidAllCondsStart[ii,:]
                                                        [listMidLightSigmoidAllConds[ii,:]==idxCommon]),ii][0][SigmoidIdx==idxCommon]
                if GaussRValue < SigmoidRValue:
                    listMidLightGaussAllConds[ii,:][listMidLightGaussAllConds[ii,:]==idxCommon] = -1
                    listMidLightGaussAllCondsStart[ii,:][listMidLightGaussAllConds[ii,:]==idxCommon] = -1
                else:
                    listMidLightSigmoidAllConds[ii,:][listMidLightSigmoidAllConds[ii,:]==idxCommon] = -1
                    listMidLightSigmoidAllCondsStart[ii,:][listMidLightSigmoidAllConds[ii,:]==idxCommon] = -1               
            sortIdxs = np.argsort(listMidLightGaussAllConds[ii,:])[::-1]            
            listMidLightGaussAllConds[ii,:] = listMidLightGaussAllConds[ii,sortIdxs]
            listMidLightGaussAllCondsStart[ii,:] = listMidLightGaussAllCondsStart[ii,sortIdxs]
            sortIdxs = np.argsort(listMidLightSigmoidAllConds[ii,:])[::-1]            
            listMidLightSigmoidAllConds[ii,:] = listMidLightSigmoidAllConds[ii,sortIdxs]
            listMidLightSigmoidAllCondsStart[ii,:] = listMidLightSigmoidAllCondsStart[ii,sortIdxs]


        cellsCommon = np.intersect1d(listHighLightSigmoidAllConds[ii,:][listHighLightSigmoidAllConds[ii,:]>-1],
                                   listHighLightGaussAllConds[ii,:][listHighLightGaussAllConds[ii,:]>-1])
        if len(cellsCommon):
            for idxCommon in cellsCommon:

                GaussIdx = dictHighLightGaussLowRSquared[int(listHighLightGaussAllCondsStart[ii,:]
                                                        [listHighLightGaussAllConds[ii,:]==idxCommon]),ii][0]
                GaussRValue = dictHighLightGaussLowRSquaredValues[int(listHighLightGaussAllCondsStart[ii,:]
                                                        [listHighLightGaussAllConds[ii,:]==idxCommon]),ii][0][GaussIdx==idxCommon]

                SigmoidIdx = dictHighLightSigmoidLowRSquared[int(listHighLightSigmoidAllCondsStart[ii,:]
                                                        [listHighLightSigmoidAllConds[ii,:]==idxCommon]),ii][0]
                SigmoidRValue = dictHighLightSigmoidLowRSquaredValues[int(listHighLightSigmoidAllCondsStart[ii,:]
                                                        [listHighLightSigmoidAllConds[ii,:]==idxCommon]),ii][0][SigmoidIdx==idxCommon]
                if GaussRValue < SigmoidRValue:
                    listHighLightGaussAllConds[ii,:][listHighLightGaussAllConds[ii,:]==idxCommon] = -1
                    listHighLightGaussAllCondsStart[ii,:][listHighLightGaussAllConds[ii,:]==idxCommon] = -1
                else:
                    listHighLightSigmoidAllConds[ii,:][listHighLightSigmoidAllConds[ii,:]==idxCommon] = -1
                    listHighLightSigmoidAllCondsStart[ii,:][listHighLightSigmoidAllConds[ii,:]==idxCommon] = -1
            sortIdxs = np.argsort(listHighLightGaussAllConds[ii,:])[::-1]            
            listHighLightGaussAllConds[ii,:] = listHighLightGaussAllConds[ii,sortIdxs]
            listHighLightGaussAllCondsStart[ii,:] = listHighLightGaussAllCondsStart[ii,sortIdxs]   
            sortIdxs = np.argsort(listHighLightSigmoidAllConds[ii,:])[::-1]            
            listHighLightSigmoidAllConds[ii,:] = listHighLightSigmoidAllConds[ii,sortIdxs]
            listHighLightSigmoidAllCondsStart[ii,:] = listHighLightSigmoidAllCondsStart[ii,sortIdxs]         
        
    """
    Plotting using average and sem the raw data and comparing against the fitted model
    Also calculating interpolated R-squared which can help us determine how close the subsampled data is to the function. 
    """

    cntGaussPositive = 0 # variable that counts how many cells fit using gauss have positive MI
    cntSigmoidPositive = 0 # variable that counts how many cells fit using sigmoid have positive MI
    x_sampled = np.arange(0,90,2)

    gaussParametersArray = np.zeros((3,1))
    for dataset in range(num_cells):
        for posCells in range(sum(listNoLightGaussAllConds[dataset,:]>-1)):
            if listNoLightGaussAllConds[dataset,posCells] > -1:
                cell = int(listNoLightGaussAllConds[dataset,posCells])
                startPt = int(listNoLightGaussAllCondsStart[dataset,posCells])
                bestParameters = dictNoLightGaussLowRSquaredFitParams[startPt,dataset][0][:,dictNoLightGaussLowRSquared
                                                                                          [startPt,dataset][0]==cell]

                gaussParametersArray = np.append(gaussParametersArray,bestParameters,axis=1)
                [amp0,mean0,std0] = bestParameters            
                avgResponse = responseFunctionMat['Ind_CellSelectionIncludingSEM'][dataset,0][[0,3,6,9,12,15,18],cell]
                semResponse = responseFunctionMat['Ind_CellSelectionIncludingSEM'][dataset,1][[0,3,6,9,12,15,18],cell]
                ifSoundPositive = check_soundPositive(parameters=bestParameters,
                                                      ifGauss=1)
                cntGaussPositive += int(ifSoundPositive)
                plt.errorbar(amplitude, avgResponse, yerr=semResponse)
                y = amp0*np.exp(-(x_sampled-mean0)**2/(2*std0**2))
                plt.plot(x_sampled, y,'k--') 
                y = amp0*np.exp(-(amplitude-mean0)**2/(2*std0**2))
                plt.plot(amplitude, y,'r--') 
                plt.xlabel('Intensity',fontsize=15)
                plt.ylabel('Response',fontsize=15)
                plt.xticks(fontsize=15)
                plt.yticks(fontsize=15)
                plt.tight_layout()
                print('params',bestParameters, 'MI',monotonicity_rawData(avgResponse), 
                      'Sound Positive?', ifSoundPositive)
                plt.show()        

    sigmoidParametersArray = np.zeros((4,1))
    interpolatedErrorArray = np.zeros((1))
    for dataset in range(num_cells):
        for posCells in range(sum(listNoLightSigmoidAllConds[dataset,:]>-1)):
            if listNoLightSigmoidAllConds[dataset,posCells] > -1:
                cell = int(listNoLightSigmoidAllConds[dataset,posCells])
                startPt = int(listNoLightSigmoidAllCondsStart[dataset,posCells])
                bestParameters = dictNoLightSigmoidLowRSquaredFitParams[startPt,dataset][0][:,dictNoLightSigmoidLowRSquared
                                                                                            [startPt,dataset][0]==cell]

                sigmoidParametersArray = np.append(sigmoidParametersArray,bestParameters,axis=1)
                [y0, ymax, x0, delta_x] = bestParameters
                avgResponse = responseFunctionMat['Ind_CellSelectionIncludingSEM'][dataset,0][[0,3,6,9,12,15,18],cell]
                semResponse = responseFunctionMat['Ind_CellSelectionIncludingSEM'][dataset,1][[0,3,6,9,12,15,18],cell]
                ifSoundPositive = check_soundPositive(parameters=bestParameters,
                                                      ifGauss=0)
                cntSigmoidPositive += int(ifSoundPositive)
                if np.abs(ymax):
                    [interpolated_err, 
                     interpolated_amp, 
                     interpolated_response] = interpolatedSigmoid_Rsquared(params=bestParameters, 
                                                                           response=avgResponse)
                    interpolatedErrorArray = np.append(interpolatedErrorArray,interpolated_err)
                    plt.errorbar(amplitude, avgResponse, yerr=semResponse)
                    plt.plot(interpolated_amp, interpolated_response,'o')
                    y = y0 + (ymax-y0)/(1+np.exp((x0-x_sampled)/delta_x))
                    valueAt20 = (y[-1]-y[0])*0.1 + y[0]
                    valueAt80 = (y[-1]-y[0])*0.9 + y[0]
                    valueAt50 = (y[-1]-y[0])*0.5 + y[0]            
                    plt.plot(x_sampled, y,'k--') 
                    y = y0 + (ymax-y0)/(1+np.exp((x0-amplitude)/delta_x))
                    plt.plot(amplitude, y,'r--') 
                    plt.plot(x0-delta_x*np.log(((ymax-y0)/(valueAt20-y0))-1),valueAt20,'r*')
                    plt.plot(x0-delta_x*np.log(((ymax-y0)/(valueAt80-y0))-1),valueAt80,'g*')
                    plt.plot(x0-delta_x*np.log(((ymax-y0)/(valueAt50-y0))-1),valueAt50,'k*')   
                    plt.xlabel('Intensity',fontsize=15)
                    plt.ylabel('Response',fontsize=15)
                    plt.xticks(fontsize=15)
                    plt.yticks(fontsize=15)
                    plt.tight_layout()
                    print('params',bestParameters, 'MI',monotonicity_rawData(avgResponse), 
                          'Sound Positive?', ifSoundPositive, 'Interpolated error', interpolated_err)
                    plt.show()        

    return [gaussParametersArray, sigmoidParametersArray, interpolatedErrorArray,
            listNoLightGaussAllConds, listMidLightGaussAllConds,
            listHighLightGaussAllConds, listNoLightSigmoidAllConds,
            listMidLightSigmoidAllConds, listHighLightSigmoidAllConds]