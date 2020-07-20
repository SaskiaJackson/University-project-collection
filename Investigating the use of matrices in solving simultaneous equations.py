#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 11:25:53 2019

@author: saskiajackson
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.linalg as scipylin

'''--------------------------TO RUN THE CODE------------------------------  '''        
    
#See below the functions for a set of commented code that can be uncommented for which sections of code want to run

'''-----------------------------PART ONE----------------------------------  '''


def randommatrix(n,m):
    return np.random.randint(9, size=(n,m))

def transpose(matrix):
    n = len(matrix)
    a = np.zeros( (n,n) )
    
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            a[j][i] = matrix[i][j]
    
    return a


def minors(matrix, i, j):
        
    return np.delete(np.delete(matrix, i, axis=0), j, axis=1)
        
def determinant(matrix, k):
    
    det = 0
    
    if len(matrix) == 2:
        return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]
        
    for k in range(len(matrix)):
        det += matrix[0][k] * determinant(minors(matrix, 0, k),k) * (-1)**(k)
    
    return det

def cofactors(matrix):
    b = np.zeros( (len(matrix),len(matrix)) )
    
    if len(matrix) == 2:
        b[0][0] = matrix[1][1]
        b[1][1] = matrix[0][0]
        b[0][1] = -matrix[0][1]
        b[1][0] = -matrix[1][0]
        return b
    
    for k in range(len(matrix)):
        for l in range(len(matrix)):
            b[k][l] = determinant(minors(matrix,k,l),k) * (-1)**(k+l)
    return b
        
def inverse(matrix):
    if len(matrix) == 2:
        return (1/determinant(matrix,0)) * cofactors(matrix)
    return (1/determinant(matrix,0)) * transpose(cofactors(matrix))

def plot(x,y, xlabel, ylabel, legend):
    plt.plot(x,y,linewidth=1, label = legend)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
def graphtime(maxn):
    sizeofmatrix = []
    timetaken = []
    for n in range(2,maxn):
        matrix = randommatrix(n,n)
        while determinant(matrix,0) == 0:
            matrix = randommatrix(n,n)
        begin = time.clock()  
        inverse(matrix)
        t = time.clock() - begin
        sizeofmatrix.append(n)
        timetaken.append(t)
        
    plot(sizeofmatrix,timetaken,'Dimensions of the matrix','Time taken to determine the inverse','Inverse of a metrix')
    plt.show()

def grapherror(maxn):
    sizeofmatrix = []
    error = []
    for n in range(2,maxn):
        matrix = randommatrix(n,n)
        while determinant(matrix,0) == 0:
            matrix = randommatrix(n,n)
        dif = np.subtract(np.matmul(inverse(matrix),matrix),np.identity(n))
        errors = np.max(dif)
        sizeofmatrix.append(n)
        error.append(errors)
    plot(sizeofmatrix,error,'Dimensions of the matrix','Maximum error of the inverse function','Inverse of a matrix')
        

'''-----------------------------PART TWO----------------------------------  '''

def solvewithinverse(matrix,b):
    return np.matmul(inverse(matrix),b)

def LUdecomposition(matrix,b):
    P,L,U = scipylin.lu(matrix)
    LU = np.matmul(np.linalg.inv(U),np.linalg.inv(L))
    LUP = np.matmul(LU,P)
    return (np.matmul(LUP,b))

def SVdecomposition(matrix,b):
    U,S,V = np.linalg.svd(matrix)
    sigma = np.diag(S)
    SU = np.matmul((np.linalg.inv(sigma)).transpose(),U.transpose())
    VSU = np.matmul(V,SU)
    return np.matmul(VSU,b)

def comparemethods(maxn):
    sizeofmatrix = []
    timetakeninverse = []
    timetakenLU = []
    timetakenSV = []
    for n in range(2,maxn):
        matrix = randommatrix(n,n)
        b = randommatrix(n,1)
        
        while determinant(matrix,0) == 0:
            matrix = randommatrix(n,n)
        
        begin = time.clock()
        solvewithinverse(matrix,b)
        inverset = time.clock() - begin
        timetakeninverse.append(inverset)
        
        begin = time.clock()  
        LUdecomposition(matrix,b)
        LUt = time.clock() - begin
        timetakenLU.append(LUt)
        
        begin = time.clock()
        SVdecomposition(matrix,b)
        SVt = time.clock() - begin
        timetakenSV.append(SVt)
        
        sizeofmatrix.append(n)
    plot(sizeofmatrix,timetakeninverse,'Dimensions of the matrix','Time taken to determine the solution','Cramer')
    plot(sizeofmatrix,timetakenLU,'Dimensions of the matrix','Time taken to determine the solution','LU decomposition')
    plot(sizeofmatrix,timetakenSV,'Dimensions of the matrix','Time taken to determine the solution','SV decomposition')
    plt.legend()
    plt.show()
    
def compareerrors(maxn):
    sizeofmatrix = []
    errorinverse = []
    errorLU = []
    errorSV = []
    for n in range(2,maxn):
        matrix = randommatrix(n,n)
        b = randommatrix(n,1)
        
        while determinant(matrix,0) == 0:
            matrix = randommatrix(n,n)
        
        difLU = np.max(np.subtract(b,np.matmul(matrix,LUdecomposition(matrix,b))))
        errorLU.append(difLU)
        
        difSV = np.max(np.subtract(b,np.matmul(matrix,SVdecomposition(matrix,b))))
        errorSV.append(difSV)
        
        sizeofmatrix.append(n)
    plot(sizeofmatrix,errorinverse,'Dimensions of the matrix','Maximum error','Cramer')
    plot(sizeofmatrix,errorLU,'Dimensions of the matrix','Maximum error','LU decomposition')
    plot(sizeofmatrix,errorSV,'Dimensions of the matrix','Maximum error','SV decomposition')
    plt.legend()
    plt.show()
    
    
def nearsingular():
    kval = []
    timetakeninverse = []
    timetakenLU = []
    timetakenSV = []
    for i in range(1,500):
        k = i/1000
        matrix = np.array([[1,1,1],
                           [1,2,-1],
                           [2,3,k]])
        b = np.array([[5],
                      [10],
                      [15]])
        begin = time.clock()
        solvewithinverse(matrix,b)
        inverset = time.clock() - begin
        timetakeninverse.append(inverset)
    
        begin = time.clock()
        LUdecomposition(matrix,b)
        LUt = time.clock() - begin
        timetakenLU.append(LUt)
        
        begin = time.clock()
        SVdecomposition(matrix,b)
        SVt = time.clock() - begin
        timetakenSV.append(SVt)
        
        kval.append(k)
    plot(kval,timetakeninverse,'Value of k in the near singular matrix','Time taken to determine the solution','Cramers Rule')
    plot(kval,timetakenLU,'Value of k in the near singular matrix','Time taken to determine the solution', 'LU decomposition')
    plot(kval,timetakenSV, 'Value of k in the near singular matrix','Time taken to determine the solution','SVdecomposition')
    plt.legend()
    plt.show()

'''----------------------------PART THREE---------------------------------  '''

def LHSanglematrix(theta1,theta2,phi,cosalpha,sinalpha,cosbeta,sinbeta,cosgamma,singamma):
    LHS = np.array([[-np.cos(theta1)*cosalpha,np.cos(theta2)*cosbeta,np.cos(phi)*cosgamma],
                     [-np.cos(theta1)*sinalpha,-np.cos(theta2)*sinbeta,np.cos(phi)*singamma],
                     [np.sin(theta1),np.sin(theta2),np.sin(phi)]])
    return LHS

def RHSanglematrix(theta1,theta2,phi,cosalpha,sinalpha,cosbeta,sinbeta,cosgamma,singamma):
    RHS = np.array([[-np.cos(theta1)*cosalpha,np.cos(theta2)*cosbeta,-np.cos(phi)*cosgamma],
                     [-np.cos(theta1)*sinalpha,-np.cos(theta2)*sinbeta,np.cos(phi)*singamma],
                     [np.sin(theta1),np.sin(theta2),np.sin(phi)]])
    return RHS
    
def tension2wires(wire):
    plot2 = [] 
    maximumT = 0
    maximumposition = []
    for Z in range(1,700):
        plot = []
        z = Z/100
        for X in range(1,1500):
            x = X/100
            theta1 = np.arctan((8-z)/(x))
            theta2 = np.arctan((8-z)/(15-x))
            thetamatrix = np.array([[np.cos(theta1),-np.cos(theta2)],
                                     [np.sin(theta1),np.sin(theta2)]])
            mass = np.array([[0],
                             [70*9.81]])
            tensionmatrix = np.matmul(np.linalg.inv(thetamatrix),mass)
            
            plot.append(float(tensionmatrix[wire]))
            if plot[-1] > maximumT:
                maximumT = plot[-1]
                maximumposition = [x,z] 
        plot2.append(plot)
    plt.imshow(plot2, origin='lower',extent=[0,15,0,7],cmap='rainbow')
    plt.xlabel('x /m')
    plt.ylabel('z /m')
    plt.colorbar(label = 'Tension /N')
    plt.show()
    print('The maximum value of tension is ',maximumT, '. This occurs at (x,z) = ',maximumposition)
                   
def tension3wires(wire):
    for z in range(1,8):
        plot2 = []
        maximumT = 0
        maximumposition = []
        for Y in range(1,800):
            plot = []
            y = Y/100
            for X in range(1,1500):
                x = X/100
                if x==7.5:
                    continue
                
                theta1 = np.arctan((8-z)/np.sqrt(x**2 + y**2))
                theta2 = np.arctan((8-z)/np.sqrt((15-x)**2 + y**2))
                phi = np.arctan((8-z)/np.sqrt((7.5-x)**2 + (8-y)**2))
                cosalpha = x/np.sqrt(x**2 + y**2)
                sinalpha = y/np.sqrt(x**2 + y**2)
                cosbeta = (15-x)/np.sqrt((15-x)**2 + y**2)
                sinbeta = y/np.sqrt((15-x)**2 + y**2)
                cosgamma = abs(x-7.5)/np.sqrt((x-7.5)**2 + (8-y)**2)
                singamma = (8-y)/np.sqrt((x-7.5)**2 + (8-y)**2)
                mass = np.array([[0],
                                 [0],
                                 [70*9.81]])
                if x < 7.5:
                    if(y >= (8/7.5)*x):
                        plot.append(0)
                    else:
                        anglematrix = LHSanglematrix(theta1,theta2,phi,cosalpha,sinalpha,cosbeta,sinbeta,cosgamma,singamma)
                        tensionmatrix = np.matmul(np.linalg.inv(anglematrix),mass)
                        plot.append(float(tensionmatrix[wire]))
                    
                elif x > 7.5:
                    if(y >= (((-8)/7.5)*x)+ 16):
                        plot.append(0)
                    else:
                        anglematrix = RHSanglematrix(theta1,theta2,phi,cosalpha,sinalpha,cosbeta,sinbeta,cosgamma,singamma)
                        tensionmatrix = np.matmul(np.linalg.inv(anglematrix),mass)
                        plot.append(float(tensionmatrix[wire]))
                if plot[-1] > maximumT:
                    maximumT = plot[-1]
                    maximumposition = [x,y]
            plot2.append(plot)
        print('-------------------  z = ', z, '  -------------------')
        
        plt.imshow(plot2, origin='lower', extent=[0,15,0,8], cmap='rainbow')
        plt.xlabel('x /m')
        plt.ylabel('y /m')
        plt.colorbar(label = 'Tension /N')
        plt.show()
        print('The maximum value of tension is ',maximumT,'. This occurs at (x,y) = ',maximumposition, 'for a height ',z)
        
'''--------------------------TO RUN THE CODE------------------------------  '''        
    
#matrix = randommatrix(3,3)    #defines matrix as a random 3x3 matrix, change the numbers in the argument to change the dimensions of the matrix 
                                           #should be run to print inverse matrix or solve simultaneous equations 
#b = randommatrix(3,1)      #defines a vector to be used in the simultaneous equations, change the numbers in the argument to change the dimensions of the matrix 
                                             #the second argument should be kept as 1 to keep it as a vector
                                             #should be run to solve simultaneous equations

#---------------------------PRINT INVERSE MATRIX-------------------------------
#print(np.linalg.inv(matrix))    #prints the inbuilt function for inverting the matrix    
#print(inverse(matrix))     #prints the made function for inverting the matrix

#-----------------------SOLVE SIMULTANEOUS EQUATIONS---------------------------
#print(solvewithinverse(matrix,b))    #prints the solution of simultaneous equations using Cramers rule to find an inverse 
#print(LUdecomposition(matrix,b))     #prints the solution of simultaneous equations using LU decomposition
#print(SVdecomposition(matrix,b))     #prints the solution of simultaneous equations using SV decomposition
 
#---------------------------GRAPHICAL SOLUTIONS--------------------------------                                           
                            #for the following three, enter in an argument to determine the maximum matrix dimension shown in the graph
graphtime(7)    #prints a graph comparing the time taken to invert the matrix using Cramers rule
grapherror(7)   #prints a graph comparing the error in the inverted matrix when Cramers rule is used                                        
#comparemethods(10)    #prints a graph comparing the time taken for the three different methods for various matrix dimensions
#compareerrors(10)

#nearsingular()       #prints a graph comparing the the time taken for the three different methods for increasingly near singular matrices
#tension2wires(0)     #prints a plot showing how the tension varies with position for two wires 
                            #the argument should be changed to determine the wire that the tension is found for; 0 for the wire from the front left drum and 1 for the wire from the front right drum
#tension3wires(2)   #prints a plot showing how the tension varies with position for three wires 
                            #the argument should be changed to determine the wire that the tension is found for; 0 for the wire from the front left drum,1 for the wire from the front right drum, 2 for the wire from the back drum


