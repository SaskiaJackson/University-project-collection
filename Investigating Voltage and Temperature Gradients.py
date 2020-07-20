#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 17:27:12 2019

@author: saskiajackson
"""
import numpy as np
import matplotlib.pyplot as plt

print(
'-------------------------------TO RUN THE CODE-------------------------------')
#below each function are a set of commented out lines that call the various functions
#to run the code, uncomment which function you would like to see


print(
'----------------------------------PART ONE-----------------------------------')


def initialconditions(n,m,left,right,bottom,top):   #n is the number of rows and m is the number of columns
    grid = 1 * np.zeros( (n,m) )   #creates grid of random integers between 0 and 50 
    for i in range(n):   
        grid[i,0] = left         #sets the values in the grid to zero along the lhs, rhs, top and bottom of the grid
        grid[i,n-1] = right
    for j in range(m):
        grid[0,j] = bottom
        grid[m-1,j] = top   
        
    return grid
    
def GaussSeidel(n,m,conditions,voltage,capacitor='0',periodicboundary='0'):
    count = 0
    sumbefore = []    
    sumafter = [] 
    grid = conditions
    percentage = 10
    while percentage > 1e-15:      #convergence condition
        sumbefore = np.copy(grid) 
        for i in range(1,n-1):       #between 1 and n-1/m-1 to cover all the points in the grid except the edges
            for j in range(1,m-1):
                if capacitor =='yes':
                    if grid[i][j] != voltage and grid[i][j] != -voltage:     #condition to ensure that the capacitor is not changed while the Gauss Seidel method is used
                        grid[i][j] = 0.25 * (grid[i-1][j] + grid[i+1][j] + grid[i][j-1] + grid[i][j+1])    
                        if periodicboundary == 'yes':
                            if i == n-1 and j != n-1:
                                grid[i][j] = (0.25)*(grid[i-1][j]+grid[0][j]+grid[i][j-1]+grid[i][j+1])
                            
                            elif j == n-1 and i != n-1:
                                grid[i][j] = (0.25)*(grid[i-1][j]+grid[i+1][j]+grid[i][j-1]+grid[i][0])
                            
                            elif i == n-1 and j == n-1:
                                grid[i][j] = (0.25)*(grid[i-1][j]+grid[0][j]+grid[i][j-1]+grid[i][0])
                            
                            else:
                                grid[i][j] = (0.25)*(grid[i-1][j]+grid[i+1][j]+grid[i][j-1]+grid[i][j+1])
                    else:
                        grid[i][j] = 0.25 * (grid[i-1][j] + grid[i+1][j] + grid[i][j-1] + grid[i][j+1])
        sumafter = np.copy(grid)
        difference = abs(sumafter - sumbefore)
        percentage = abs(np.amax(difference)/np.amax(sumafter))      #defines a percentage difference of the matrix before and after iteration 
        count += 1     #determines the number of iterations completed until the convergence condition is satisfied
    print('Number of iterations = ', count)
    return grid

def colourplot(n,m,array,Xlabel,Ylabel,title,barlabel,tmax,arrows='0',extent='0'):
    X,Y = np.meshgrid(np.arange(0,n), np.arange(0,m))    #creates a mesh grid of points that can be plotted onto
    if extent == 'yes':
        plt.contourf(X,Y, array,100,cmap='jet',extent=[0,tmax,0,0.5])
    else:
        plt.contourf(X,Y, array,100,cmap='jet')      #plots the mesh grid with the grid values in a contour plot
    plt.colorbar(label=barlabel)
    plt.title(title)
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)
    
    if arrows == 'yes':      #adds arrows on the colour plot to indicate the electric field
        d = np.gradient(array)
        plt.streamplot(X, Y, -d[1], -d[0], density = 1, color = 'black')
    plt.show()
    
def capacitorboundary(n,m,sep,length,voltage):
    grid = np.zeros( (n,m) )
    top = int((m + sep)/2)        #position of the top capacitor
    left = int((n - length)/2)     #position of the left hand side of the capacitors
    right = int((n + length)/2)     #position of the right hand side of the capacitors
    bottom = int((m - sep)/2)        #position of the bottom capacitor 
    grid[top,left:right] = voltage       #sets the top capacitor to 100
    grid[bottom,left:right] = -voltage
    return grid

#def AdjacentPoints(V,i,j,Periodic=False):
#    #This function averages over the adjacent 4 point in the grid 
#    #If periodic is true at the edge of the grid it will take the value of the 
#    #grid at the opposite side.
#    n=V.shape
#    if Periodic==False:
#        #There is a subtelty that I had to keep in mind when progrmaing this
#        #When i or j is equal to 0, i-1 or j-1 will be -1 which will return the 
#        #opposite bouandry point without any additional if statments.
#        if i==0 or j==0 or i==n[0]-1 or j==n[1]-1:
#            return V[i,j]
#        else:
#            return (1/4)*(V[i-1,j]+V[i+1,j]+V[i,j-1]+V[i,j+1])
#    elif Periodic==True: 
#        if i==n[0]-1 and j!=n[1]-1:
#            return (1/4)*(V[i-1,j]+V[0,j]+V[i,j-1]+V[i,j+1])
#        
#        elif j==n[1]-1 and i!=n[0]-1:
#            return  (1/4)*(V[i-1,j]+V[i+1,j]+V[i,j-1]+V[i,0])
#        
#        elif i==n[0]-1 and j==n[1]-1:
#            return (1/4)*(V[i-1,j]+V[0,j]+V[i,j-1]+V[i,0])
#        
#        else:
#            return (1/4)*(V[i-1,j]+V[i+1,j]+V[i,j-1]+V[i,j+1])
    
def compare():
     xaxis = np.linspace(0,100,20)
     x = [10,20,30,40,50,60,70,80,90,100]
     y = [257,1099,2490,4415,6855,9808,13267,17216,21662,26591]    
     z = np.polyfit(x,y,2)
     p = np.poly1d(z)
     plt.plot(x,y,'.',xaxis,p(xaxis),'--')
     plt.xlabel('Number of points on the grid')
     plt.ylabel('Number of iterations')
     plt.title('Comparing how the number of iterations change as the grid size increases')
     plt.show()
     
n = 20
m = n
voltage = 50
tmax=0

#colourplot(n,m,GaussSeidel(n,m,initialconditions(n,m,0,0,0,100),voltage),' ',' ','Voltage over a known area with a set value of 100V along the top edge.','Voltage /V',tmax)
#
#colourplot(n,m,GaussSeidel(n,m,capacitorboundary(n,m,4,10,voltage),voltage,capacitor='yes'),' ',' ','Voltage experienced from two capacitor plates with opposite voltage.','Voltage /V',tmax,arrows='yes',)

colourplot(n,m,GaussSeidel(n,m,capacitorboundary(n,m,4,n,voltage),voltage,capacitor='yes',periodicboundary='yes'),' ',' ','Voltage experienced from two capacitor plates with opposite voltage.','Voltage /V',tmax,arrows='yes',)

#compare()

print(
'----------------------------------PART TWO-----------------------------------')

def definematrix(h, dt ,n, ice='0'):      #h is the space step, dt is the time step, n is the number of points on the array
    alpha = 59/(450*7900)       #alpha = thermal conductivity/(density * specific heat)
    x = (alpha*dt) / (h**2)     #defines a value for x to use in the matrix
    matrix = np.zeros( (n,n) )
    for i in range(1,n-1):
        matrix[i][i] = 1+2*x      #sets the diagonal elements of the matrix to 1 + 2x
        matrix[i][i+1] = -x       #sets the off diagonal elements of the matrix to -x
        matrix[i][i-1] = -x
    
    matrix[0][0] = 1+3*x          #sets the top left value of the matrix to be 1 + 3x
    matrix[-1][-1] = 1+x          #sets the bottom right corner of the matrix to be 1 + x
    matrix[0][1] = matrix[-1,-2]=-x      #sets the final off diagonal elements to -x
    
    if ice == 'yes':        #used when one end of the rod is emersed in an ice bath at 0 degrees
        matrix[-1][-1] = 1 + 3*x      #sets the bottom right corner of the matrix to be 1 + 3x (the same as the top left corner as they are now both boundary conditions)
    
    return matrix,x
        
def temparray(n):     #creates an array to for the inial rod matrix
    rod = 20*np.ones( (n,1) )
    return rod

def diffusion(n, tmax,Nt,ice='0'):      #tmax is the maximum time that the simulation is run for
    dt = tmax/Nt        #determines the timestep 
    h=0.5/n        #determines the value of h from the length of the rod (50cm) and the number of point in the array 
    matrix,x = definematrix(h, dt, n)
    if ice == 'yes':
        matrix,x = definematrix(h, dt, n, ice='yes')
    rod = temparray(n)
    rodarray = np.linspace(0,0.5,n)     #used in plotting the length of the rod
    
    for t in range(Nt):     
        rod[0][0] += 2*x*1000       #sets the initial value of the rod to 2*x*voltage of the heat reservoir
        rod = np.linalg.solve(matrix,rod)     #solves the matrix equation to find the rod temperatures after one iteration 
        if t%(Nt/10) == 0:
            lineplot(rodarray,rod,int(tmax*(t/Nt)),'Distance along the rod /m','Temperature /degrees celsius','A graph comparing the temperature along the rod at various times')
    plt.show()
    return rod
    
def comparisongraph(n,tmax,Nt,ice='0'):      #tmax is the maximum time that the simulation is run for
    dt = tmax/Nt        #determines the timestep 
    h=0.5/n        #determines the value of h from the length of the rod (50cm) and the number of point in the array 
    matrix,x = definematrix(h, dt, n)
    if ice == 'yes':
        matrix,x = definematrix(h, dt, n, ice='yes')
    rod = temparray(n)
    comparison = rod.transpose()
    
    for t in range(Nt):
        rod[0][0] = 2*x*1000
        rod = np.linalg.solve(matrix,rod)     #solves the matrix equation to find the rod temperatures after one iteration 
        b = rod.transpose()
        comparison = np.concatenate((comparison,b))
    comparison = np.delete(comparison,0,axis=0)
    colourplot(n,n,comparison,'Distance along rod','Time taken','Comparison of how the temperature changes along the rod and the time increases','Temperature /degrees celcius',tmax,extent='yes')
    
def lineplot(Xarray,Yarray,label,Xlabel,Ylabel,title):
    plt.plot(Xarray,Yarray, label=label)
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)
    plt.title(title)
    plt.legend()
    
def gradientplot(array,Ylabel,title,bartitle):
    plt.cur_axes = plt.gca()
    plt.cur_axes.axes.get_xaxis().set_ticks([])
    plt.imshow(array, cmap='jet', vmin=0, extent = [0,17,50,0])
    plt.ylabel(Ylabel)
    plt.title(title)
    plt.colorbar(label=bartitle)
    plt.show()
    
    
n=300           #number of points in space used
tmax = 20000    #maximum time the simulation is run for 
tmax2 = 3600    #the value of tmax used in the snapshot of temperature gradient over the rod distance
Nt = 1000       #number of points in time used

#gradientplot(diffusion(n,tmax2,Nt),'Distance along rod /cm','How the temperature varies along the rod at time %is' %tmax2,'Temperature /degrees celcius')
#print('When there is no heat loss from the far end of the rod')
#diffusion(n,tmax,Nt)
#comparisongraph(n,tmax,300)
#print('When one side of the rod is imersed in a block of ice at 0 degrees celcius')
#diffusion(n,tmax,Nt,ice='yes')
#comparisongraph(n,tmax,300,ice='yes')


