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


amplitude = np.array([0,30,50,60,70,80,90],dtype=np.float128)

def load_ResponseData(responses_expt,sem_expt):
    response_noLightOrig = np.array(responses_expt[[0,3,6,9,12,15,18],:],dtype=np.float128)
    response_midLightOrig = np.array(responses_expt[[1,4,7,10,13,16,19],:],dtype=np.float128)
    response_highLightOrig = np.array(responses_expt[[2,5,8,11,14,17,20],:],dtype=np.float128)    
    response_noLightNorm = np.zeros_like(response_noLightOrig)
    response_midLightNorm = np.zeros_like(response_midLightOrig)
    response_highLightNorm = np.zeros_like(response_highLightOrig)
    
    sem_noLightOrig = np.array(sem_expt[[0,3,6,9,12,15,18],:],dtype=np.float128)
    sem_midLightOrig = np.array(sem_expt[[1,4,7,10,13,16,19],:],dtype=np.float128)
    sem_highLightOrig = np.array(sem_expt[[2,5,8,11,14,17,20],:],dtype=np.float128)
    sem_noLightNorm = np.zeros_like(sem_noLightOrig)
    sem_midLightNorm = np.zeros_like(sem_midLightOrig)
    sem_highLightNorm = np.zeros_like(sem_highLightOrig)
    
    numSoundResponsiveCells = responses_expt.shape[1]
    for iCell in range(numSoundResponsiveCells):
        response_noLightNorm[:,iCell] = response_noLightOrig[:,iCell]/np.max(np.abs(response_noLightOrig[:,iCell]))
        response_midLightNorm[:,iCell] = response_midLightOrig[:,iCell]/np.max(np.abs(response_midLightOrig[:,iCell]))  
        response_highLightNorm[:,iCell] = response_highLightOrig[:,iCell]/np.max(np.abs(response_highLightOrig[:,iCell])) 
        sem_noLightNorm[:,iCell] = sem_noLightOrig[:,iCell]/np.max(np.abs(response_noLightOrig[:,iCell]))
        sem_midLightNorm[:,iCell] = sem_midLightOrig[:,iCell]/np.max(np.abs(response_midLightOrig[:,iCell]))
        sem_highLightNorm[:,iCell] = sem_highLightOrig[:,iCell]/np.max(np.abs(response_highLightOrig[:,iCell]))
        
    return [response_noLightNorm, response_midLightNorm, response_highLightNorm, 
            sem_noLightNorm, sem_midLightNorm, sem_highLightNorm]

def gaussian_fn(params,response,sem):
    """
    Error is based on deviance.
    """
    amp0, mean0, std0, shift = params[0], params[1], params[2], params[3]
    y = shift + amp0*np.exp(-(amplitude-mean0)**2/(2*std0**2))
    err = np.sum((y-response)**2/(2*sem**2))/np.sum((response-np.mean(response))**2/(2*sem**2))
    #err = np.sum((y-response)**2)/np.sum((response-np.mean(response))**2)
    return err

def sigmoid_fn(params,response,sem):
    """
    Error is based on deviance.
    """
    y0, ymax, x0, delta_x = params[0], params[1], params[2], params[3]
    y = y0 + (ymax-y0)/(1+np.exp((x0-amplitude)/delta_x));
    err = np.sum((y-response)**2/(2*sem**2))/np.sum((response-np.mean(response))**2/(2*sem**2))
    #err = np.sum((y-response)**2)/np.sum((response-np.mean(response))**2)
    return err

def monotonicity_gaussianFit(params,response):
    """
    Using the monotonicity index of paper - Level-Tuned Neurons in Primary Auditory Cortex Adapt Differently to Loud versus Soft
Sounds
    """
    amp0, mean0, std0, shift = params[0], params[1], params[2], params[3]
    y = shift + amp0*np.exp(-(amplitude-mean0)**2/(2*std0**2))
    mi = np.abs(y[-1]-y[0])/np.abs(np.max(y)-np.min(y))
    return mi

def monotonicity_rawData(response):
    """
    Using the monotonicity index of paper - Level-Tuned Neurons in Primary Auditory Cortex Adapt Differently to Loud versus Soft
Sounds
    """
    mi = np.abs(response[-1]-response[0])/np.abs(np.max(response)-np.min(response))
    return mi

def check_soundPositive(parameters, ifGauss):
    """
    deciding on sound positive cells - if gauss they would have positive amplitude, and if sigmoid they would have 
    either ymax>y0 and positive slope or ymax<y0 and negative slope
    """
    if ifGauss:
        ifSoundPositive = parameters[0]>0
    else:
        ifSoundPositive = (parameters[3]>0 and parameters[1] > parameters[0]) or (parameters[3]<0 and parameters[0] > parameters[1])
        if ifSoundPositive and parameters[3]<0:
            swap = np.copy(parameters[0])
            parameters[0] = np.copy(parameters[1])
            parameters[1] = swap
            parameters[3] = -parameters[3]
        elif ~ifSoundPositive and parameters[3]>0:
            swap = np.copy(parameters[0])
            parameters[0] = np.copy(parameters[1])
            parameters[1] = swap
            parameters[3] = -parameters[3]  
    return ifSoundPositive, parameters

def interpolatedSigmoid_Rsquared(params, response):
    """
    deciding what would be the cutoff and the bounded optimization method to be used for sigmoid fits. As of 06/04/21 gauss fits are yet to be tested. 
    """
    y0, ymax, x0, delta_x = params[0], params[1], params[2], params[3]
    yfit = y0 + (ymax-y0)/(1+np.exp((x0-amplitude)/delta_x));
    interpolated_amp = np.zeros((len(amplitude)*2-1,))
    interpolated_yfit = np.zeros((len(amplitude)*2-1,))
    for ii in range(len(interpolated_amp)):
        if ii%2==0:
            interpolated_amp[ii] = amplitude[ii//2]
            interpolated_yfit[ii] = yfit[ii//2]
        else:
            interpolated_amp[ii] = np.mean(amplitude[ii//2:ii//2+2])
            interpolated_yfit[ii] = np.mean(yfit[ii//2:ii//2+2])
    y = y0 + (ymax-y0)/(1+np.exp((x0-interpolated_amp)/delta_x));
    err = np.sum((y-interpolated_yfit)**2)/np.sum((interpolated_yfit-np.mean(interpolated_yfit))**2)
    return err, interpolated_amp, interpolated_yfit

def interpolatedGauss_Rsquared(params, response):
    """
    deciding what would be the cutoff and the bounded optimization method to be used for sigmoid fits. As of 06/04/21 gauss fits are yet to be tested. 
    """
    amp0, mean0, std0, shift = params[0], params[1], params[2], params[3]
    yfit = shift + amp0*np.exp(-(amplitude-mean0)**2/(2*std0**2))
    interpolated_amp = np.zeros((len(amplitude)*2-1,))
    interpolated_yfit = np.zeros((len(amplitude)*2-1,))
    for ii in range(len(interpolated_amp)):
        if ii%2==0:
            interpolated_amp[ii] = amplitude[ii//2]
            interpolated_yfit[ii] = yfit[ii//2]
        else:
            interpolated_amp[ii] = np.mean(amplitude[ii//2:ii//2+2])
            interpolated_yfit[ii] = np.mean(yfit[ii//2:ii//2+2])
    y = shift + amp0*np.exp(-(interpolated_amp-mean0)**2/(2*std0**2))
    err = np.sum((y-interpolated_yfit)**2)/np.sum((interpolated_yfit-np.mean(interpolated_yfit))**2)
    return err, interpolated_amp, interpolated_yfit
    
    
    