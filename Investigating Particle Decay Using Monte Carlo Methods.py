#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 13:55:38 2019

@author: saskiajackson
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm
import sys
import time
from scipy.stats import chisquare

#To run the code:
                 #the functions, split into sections for each part of the exercise, are at the top
                 #then all the functions are called in order to print every graph after the functions
                 

'''-----------------------------FUNCTIONS---------------------------------  '''

'''------------------------PLOTTING FUNCTIONS-----------------------------  '''

def plot(xarray,yarray,title,xlabel,ylabel,legend,scatter='0'):     #simple plot function
    if scatter == 'yes':      #optional function to plot a scatter plot instead of a line graph
        plt.plot(xarray,yarray,'.', label=legend)
    else:
        plt.plot(xarray,yarray, label=legend)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
def histogramplot(array,title,xlabel,ylabel):        #simple one dimensional histogram
    plt.hist(array, bins='auto', normed=True)        #bins set to auto so that the computer figures out the best bin size for the data provided
    plt.title(title)                                 #normed=True meaning that the axes of the histogram are normalised to enable comparison
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
def twodhistogram(xarray,yarray,title):         #2D histogram plot function
    plt.hist2d(xarray,yarray, bins = 50, norm=LogNorm())      #bins set to 50    
    plt.axis('scaled')                                        #scaled axis 
    plt.colorbar(label='Number of particles')                 #label the colour bar 
    plt.title(title)
    plt.xlabel('x position on the detector /m')
    plt.ylabel('y position on the detector /m')
    plt.show()
    
'''-----------------------------PART ONE----------------------------------  '''

def analytical(iteration):     #function for finding the distribution of the sine curve using the analytical method
        
    axis = []
    
    for i in range(iteration):
        a = np.random.uniform(0,1)     #creates a random number between 0 and 2
        axis.append(np.arccos(1-2*a))
        
    return axis          #returns the array that is plotted in the histogram
    
def montecarlo(iteration):       #function for finding the distribution of the sine curve using the accept-reject method
    
    axis = []
    
    for i in range(iteration):
        x = np.random.uniform(0,np.pi)
        y = np.random.random()
        if y < np.sin(x):         #checks if y is less than sin(x) in which case is appended to the array to be plotted in the histogram
            axis.append(x)
        
    return axis
    
def plotsin(iteration):       #function that plots the two approximations for the sine curve 
    
    X = np.arange(0,3.14,0.1)       
    Y = 0.5 * np.sin(X)           #used to plot the numpy sine function to compare the approximation to
    
    analyticalaxis = analytical(iteration)
    montecarloaxis = montecarlo(iteration)

    histogramplot(analyticalaxis,'A histogram showing the analytical method of determining a sine curve','Angle /radians','Frequency')
    plt.plot(X,Y)
    plt.show()

    histogramplot(montecarloaxis,'A histogram showing the accept-reject method of determining a sine curve','Angle /radians','Frequency')
    plt.plot(X,Y)
    plt.show()
    
def comparisontime(maxit,meaniteration):          #function to compare the time taken for each method against the number of iterations for both methods
    
    itaxis = []        #array for the number of iterations to be appended to
    Ataxis = []        #array for the time taken for the analytical/accept-reject method to be appended to
    MCtaxis = []       
    Atestlist = []     #an array to call the analytical/accept-reject function to
    MCtestlist = []    
    
    print('The graph to compare time taken with the number of iterations takes some time to compute. Please wait until the timer gets to ', maxit)
    
    for it in range(0,maxit,int(maxit/10)):
        
        sys.stdout.write('\r' + str(it))      #sets a timer to print out the number of iterations
       
        for j in range(meaniteration):
            A=0      #resetting the original time taken as 0
            MC=0
            
            begin = time.perf_counter()
            Atestlist.append(analytical(it))
            At = time.perf_counter() - begin             #measures the time that it takes to run the method
            A += At                                      #adds up the times so that the mean can be taken 
            
            start = time.perf_counter()
            MCtestlist.append(montecarlo(it))
            MCt = time.perf_counter() - start
            MC += MCt
        
        Atmean = A/meaniteration           #calculating the mean
        MCtmean = MC/meaniteration
    
        itaxis.append(it)
        Ataxis.append(Atmean)
        MCtaxis.append(MCtmean)
        
    plot(itaxis,Ataxis, 'A graph showing how the time taken for each method varies as the number of iterations is varied', 'Number of iterations', 'Time /s','Analytical method',scatter='yes')
    plot(itaxis,MCtaxis, 'A graph showing how the time taken for each method varies as the number of iterations is varied', 'Number of iterations', 'Time /s','Accept-Reject method',scatter='yes')
    plt.legend()
    plt.show()

def chisquared(method,iteration):           #function to determine the chi squared value
    
    axis = method     #method as an argument of the function so that the chi squared function can be used for both methods
    
    p,q = np.histogram(axis,100, (0,np.pi), density=True)       #p=values of the histogram on yaxis, q=edges of each bin on the xaxis
    y = 0.5 * np.sin(q[0:100])       #y is the y axis values of the numpy sine function      #q[0:100] is used to ensure that the sizes of the arrays are the same due to the extra value in the p array as the edge of the bin is used
    
    a,b = chisquare(p[1:-1],y[1:-1])       #in built function to determine the chi squared value    #a=chi squared value, b= p-value (not used)
    return a

'''-----------------------------PART TWO----------------------------------  '''

def gamma(numberparticles,r):           #function to determine where the decay occurs and to show the decay being isotropic
    
    axis = []   #axis used in the 1D histogram 
    dx = []      #x,y&z axes used in the 3D plot
    dy = []
    dz = []
    vel = 2000                 #velocity of the nuclei
    tau = 550e-6               #lifetime of the nuclei
    
    for i in range(numberparticles):
        dist = np.random.uniform(0,2)        #distance from the injection point at which the particle decays is randomly generated
        N = np.exp(-dist/(vel*tau))          #nuclear decay equation to find the number of particles decayed
        y = np.random.random()

        theta = np.arccos(1 - 2*np.random.random())       #theta and phi are the angles that determine the movement of the gamma ray after the decay
        phi = np.random.uniform(0,2*np.pi)
        
        if y < N:
            axis.append(dist)
            
            dxVal = r * np.sin(theta) * np.cos(phi)          #position of the gamma ray in sperical coordinates
            dyVal = r * np.sin(theta) * np.sin(phi)
            dzVal = r * np.cos(theta)
            dx.append(dxVal)
            dy.append(dyVal)
            dz.append(dzVal)
    
    histogramplot(axis,'A Histogram showing the distribution along the z axis of the position of the nuclei decay','Distance from the injection point /m','Frequency of nuclei decaying')        
    plt.show()
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')         #allows you to plot a 3D plot
    ax.scatter(dx,dy,dz, s=0.01)
    ax.set_aspect('equal')           #set each axis to be equal
    plt.title('A 3D plot showing the isotropic distributions of the gamma radiation')
    ax.set_xlabel('x /m')
    ax.set_ylabel('y /m')
    ax.set_zlabel('z /m')
    plt.show()
    
    
def detector(numberparticle,detectorsize,sd):               #function to determine the positions of the gamma particles on a detector screen
    
    xaxis = []
    yaxis = []
    vel = 2000                 #velocity of the nuclei
    tau = 550e-6               #lifetime of the nuclei
    
    for i in range(numberparticle):
        dist = np.random.uniform(0,2)          #distance from the injection point at which the particle decays is randomly generated
        N = np.exp(-dist/(vel*tau))            #nuclear decay equation to find the number of particles decayed
        y = np.random.random()
        r = 2 - dist                           #distance from the detector to the position of decay
        
        theta = np.random.uniform(0,2*np.pi)        #angles between the position of decay and position on the detector screen
        phi = np.arccos(1 - 2*np.random.random())
        
        if y < N and abs(r*np.cos(theta)/np.tan(phi)) < detectorsize and abs(r*np.sin(theta)/np.tan(phi)) < detectorsize and theta < np.pi:          #conditions: if position of decay appears on the detector screen in both the x and y direction
            xaxis.append(r*np.cos(theta)/np.tan(phi))          #appending the position of the decay to the arrays
            yaxis.append(r*np.sin(theta)/np.tan(phi))
            
    twodhistogram(xaxis,yaxis,'A 2D histogram showing the intensity of particles at each point on the detector')
    plt.show()
    
    xaxisSmear = np.random.normal(loc=xaxis, scale = 0.1/sd)  #0.1 = 10cm = x resolution of detector    #applies the smearing to the distribution of points on the detector 
    yaxisSmear = np.random.normal(loc=yaxis, scale = 0.3/sd)  #0.3 = 30cm = y resolution of detector  
   
    twodhistogram(xaxisSmear,yaxisSmear,'A 2D histogram showing the intensity of particles at each point on the detector once smearing has been applied')
    plt.show()

'''----------------------------PART THREE---------------------------------  '''    

def confidencelevels(iterations,maxsigma):           
    
    xaxis = []
    yaxis = []
    hline = []
    L = 12                     #luminosity
    
    print('The graph to calculate the percentage confidence for different cross sections takes some time to compute. Please wait until the timer gets to ', maxsigma)
    
    for Sigma in np.arange(0,maxsigma,0.01):       #varying cross section between o and max 
        countevents = 0
        
        sys.stdout.write('\r' + str(Sigma))    #sets a timer to print out sigma values 
        
        for i in range(iterations):
            lambdaBg = np.random.normal(5.7,0.4)       #creates a random Gaussian distribution around 5.7+-0.4
            lambdaSignal = L * Sigma                   #expectation value of the signal
            NBg = np.random.poisson(lambdaBg)          #creates a poisson distribution around the expectation value 
            NSignal = np.random.poisson(lambdaSignal)
            NEvents = NBg + NSignal
            
            if NEvents > 5:          
                countevents += 1
            
            percentage = (countevents / iterations) * 100    #number of significant events divided by the total number of events to calculate the percentage confidence 
            
        xaxis.append(Sigma)
        yaxis.append(percentage)
        hline.append(95)
    
    for i in range(0,len(yaxis)):         #used to determine the cross section value that corresponds to 95% confidence 
        if yaxis[i] > 95:
            print('')
            print('The cross section at which the confidence level is at 95% is ', xaxis[i],'nb')
            break
    
    plot(xaxis,yaxis,'A graph of the percentage confidence of the measurement for different cross sections','Cross section /nb','Percentage confidence /%','Percentage confidence', scatter='yes')
    plot(xaxis,hline,'A graph of the percentage confidence of the measurement for different cross sections','Cross section /nb','Percentage confidence /%','95% confidence')
    plt.legend()
    plt.show() 

            
'''-------------------------RUNNING THE CODE------------------------------  '''    

print(
'----------------------------------PART ONE-----------------------------------')
iteration = 100000      #number of iterations of the method
plotsin(iteration)      #input number of iterations 
plt.show()
comparisontime(100000,10)      #input the maximum number of iterations for which you want to investigate the time taken for and the number of repeats a mean is taken over
plt.show()
print('The chi squared value from the analytical method: ',chisquared(analytical(iteration),iteration))
print('The chi squared value from the accept-reject method: ',chisquared(montecarlo(iteration),iteration))
print(
'----------------------------------PART TWO-----------------------------------')
gamma(100000,1)        #input number of nuclei and radius of the volume in which the gamma particles are detected to determine if isotropic
plt.show()
detector(100000,10,1)    #input number of nuclei, detector size and standard deviation
plt.show()
print(
'---------------------------------PART THREE----------------------------------')
confidencelevels(100000,0.6)
plt.show()