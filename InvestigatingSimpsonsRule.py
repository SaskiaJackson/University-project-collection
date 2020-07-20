#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 17:37:32 2017

@author: saskiajackson
"""

import math as mt
import numpy as np
import matplotlib.pyplot as plt
import cmath 

def simpsonsrule (func, a, b, n): 
    c = (b - a)/n #the width of the function split into segments
    p = func(a) + func(b)
    
    for t in range(1, n, 2):  #sum of the odd term functions
        p+= 4*func((t*c)+a)
    for t in range(2, (n-1), 2): #sum of the even term functions
        p+= 2*func((t*c)+a)
    
    return p*(c/3)   #returns the full sum multiplied by the height


def func(X,x):   #defines a function for X(x,y',z), inside the integral of equation 5 in the questions
    r = ((1j*k/(2*z))*((X-x)**2))
    s = cmath.exp(r)   #uses the cmath module due to the complex number in the exponential
    return s


def simpsons1d(func,a,b,X,n):   #defines a new simpsons rule function but with x as a variable in order to allow x to be changed
    dx = (b - a)/n
    p = func(a, X) + func(b,X)
    
    for t in range(1, n, 2):  #sum of the odd term functions 
        p+= 4*func((t*dx)+a, X)
    for t in range(2, (n-1), 2):   #sum of the even term functions 
        p+= 2*func((t*dx)+a, X)
        
    return p*(dx/3)


def XElectric(X,x):  #a function for the inside of equation 5
    r = np.exp((1j*k/(2*z))*((X-x)**2)) 
    return r

def Electric(Y, y):   #a function describing the inside of the integral for electric field 
    R = ((1j*k/(2*z))*((Y-y)**2))
    return (simpsons1d(XElectric,a,b,X,n)*(np.exp(R)))

def XElectricC(x):   #a function for the circle aperture
    r = np.exp((1j*k/(2*z))*((X-x)**2)) 
    return r

def ElectricT(X,Y):    #the function to describe the electric field for the triangular aperture
    m = 0.0001   
    R = ((1j*k/(2*z))*((Y-y)**2))
    return simpsons1d(XElectric,-m,m,X,n)*(np.exp(R))   #integrating X(x,y,z) then multiplying by E^(ik/2z*(x.x')^2)
   

MyInput ='0'
while MyInput != 'q':
    MyInput = input('''Enter a choice:  "a" to use Simpsons rule to evaluate the integral of sin(x),
                 "b" for Fresnel diffraction in one dimension,
                 "c" for Fresnel diffractions in two dimensions from a square aperture or
                 "q" to quit: ''')  #menu system
    print('You entered the choice: ', MyInput)
    if MyInput == 'a':
        print('You have chosen part (a)')
        
        MyInput1 = input("Enter an even value for N, the number of iterations of simpson's rule: ",) 
                   #in order for the user to input a value for N
            
        try:
            n = int(MyInput1)  #checking the inputted value for N is an integer
        except ValueError:
            print("The value you entered for N is not an integer, please try again.")
            
        else:
            if n%2 == 0 :  #checking the input for N is even
                a = 0   #defining values for a, b and n
                b = mt.pi
                n = int(MyInput1)
                
                print("The integral of sin(x) using simpsons rule with", n, "iterations is ",simpsonsrule(mt.sin, a, b, n))
            else:
                print("The value you entered for N is not even, please try again.")
        
    elif MyInput == 'b':
        print('You have chosen part b')
    
        MyInput2 = input("Enter a value for N, the number of iterations of the functions: " )
        MyInput3 = input("Enter a value for z, the distance between the aperture and the screen, in millimetres: ")
                    #user imputted values for N and z
        
        try:
            n = int(MyInput2)   #to check the input for N is an integer
        except ValueError:
            print("The value you entered for N is not an integer, please try again.")
        else:
            try:
                zm = int(MyInput3)    #to check the input for z is an integer
            except ValueError:
                print("The value you entered for z is not an integer, please try again.")
            else:
                a = -0.0001  #sets aperture size
                b = 0.0001
                w = 200E-9      #wavelength
                k = 2*np.pi/w   #wavenumber
                z = zm * 1E-3   #allows the inputted value to be in mm
                
                NpPoints = 100
                x1 = -0.0001    #screen size
                x2 = 0.0001
                xaxis = np.linspace(x1,x2,NpPoints)
                yaxis = np.zeros(NpPoints)
        
                for i in range(0,NpPoints):
                    yaxis[i] = (np.absolute(simpsons1d(func,a,b,xaxis[i],n)))**2
                            #xaxis[i] replaces x in order to plot the x values on the x axis
        
           
                plt.plot(xaxis,yaxis)  #plots the graph
                plt.xlabel("x")
                plt.ylabel("|X(x,y',z)|^2")
                plt.show()

        
        
        
    elif MyInput == 'c':
        print('You have chosen part c')
        
        MyInput4 = input("Enter a value for z, the distance between the aperture and the screen, in millimetres: ")
    
        try:
            zm = float(MyInput4)  
        except ValueError:
            print("The value you entered for z is not an integer, please try again.")
        else:
            a = -0.00009    #aperture size
            b = 0.00009
            m = 0.00005     #screen size
            x1 = -m
            x2 = m
            y1 = -m
            y2 = m
            
            
            w = 400E-9         #wavelength
            k = (2*np.pi/w)    #wavenumber
            z = zm * 1E-3      #allows the input for z to be in mm
            n = 50
            E0 = 1             
            epsilon = 8.85418782E-12     #permittivity of free space
            c = 3E8                      #speed of light

 
            NpPoints = 100
            dx = (x2-x1)/n   
            dy = (y2-y1)/n
            intensity = np.zeros( (NpPoints, NpPoints) )  #sets up an inital blanck array for the intensity graph
            
            for i in range(NpPoints):
                X = x1*2 + (i * dx)       #sets the code to run through all the points on the x axis
                for j in range(NpPoints):
                    Y = y1*2 + (j * dy)   #sets the code to run through all the points on the y axis
                    Elec=simpsons1d(Electric,a,b,Y,n)   #call function here to prevent having to call twice in intensity equation
                    intensity[i,j] = epsilon * c * (abs(Elec))**2   #equation to determine the intesity 
            plt.imshow(intensity)    #plots intesity graph
            plt.colorbar()           #creates a scale next to the graph
            plt.show()
        
        
    elif MyInput == 'd':
        print('You have chosen part d')
        
    
        MyInput5 = input('Enter a choice; "c" for circle, "t" for triangle and "e" to return to the main option screen:')
        if MyInput5 == 'c':
            
            MyInput4 = input("Enter a value for z, the distance between the aperture and the screen, in millimetres: ")
               
            try:
                zm = float(MyInput4)
            except ValueError:
                print("The value you entered for z is not an integer, please try again.")
            else:
                y = float()
                rad = 0.0001    #radius of the circle
                circ = cmath.sqrt(rad**2 - y**2)   #equation of the circle
                m = 0.0001    #screen size
                x1 = -m
                x2 = m
                y1 = -m
                y2 = m
            
                w = 400E-9     
                k = (2*np.pi/w)
                z = zm * 1E-3
                n = 50
                E0 = 1
                epsilon = 8.85418782E-12
                c = 3E8
                
                def ElectricC(y):    #function for circular aperture
                    circ = np.sqrt(rad**2 - y**2)
                    R = ((1j*k/(2*z))*((Y-y)**2))
                    return simpsonsrule(XElectricC,-circ,circ,n)*(np.exp(R))
                
                NpPoints = 100
                dx = (2*m)/n
                dy = (2*m)/n
                intensity = np.zeros( (NpPoints, NpPoints) )  #sets up an inital blanck array for the intensity graph
                
                for i in range(NpPoints):
                    X = x1*2 + (i * dx)
                    for j in range(NpPoints):
                        Y = y1*2 + (j * dy)
                        Elec=simpsonsrule(ElectricC,-circ,circ,n)
                        intensity[i,j] = np.real(epsilon * c * (abs(Elec))**2)
                plt.imshow(intensity)
                plt.colorbar()
                plt.show()
            
           
        elif MyInput5 == 't':
            
            MyInput4 = input("Enter a value for z, the distance between the aperture and the screen, in millimetres: ")
    
            try:
                zm = float(MyInput4)
            except ValueError:
                print("The value you entered for z is not an integer, please try again.")
            else:
                m = 0.0001    #screen size
                x1 = -m
                x2 = m
                y1 = -m
                y2 = m
            
            
                w = 400E-9          #wavelength
                k = (2*np.pi/w)     #wavenumber
                z = zm * 1E-3      #allows the input for z to be in mm
                n = 50
                E0 = 1
                epsilon = 8.85418782E-12      #permittivity of free space 
                c = 3E8                       #speed of light

 
                NpPoints = 100
                dx = (x2-x1)/n
                dy = (y2-y1)/n
                intensity = np.zeros( (NpPoints, NpPoints) )   #sets up an inital blanck array for the intensity graph
                
                for i in range(NpPoints):    #sets the code to run through all the points on the x axis
                    X = x1*2 + (i * dx)
                    for j in range(NpPoints):   #sets the code to run through all the points on the y axis
                        Y = y1*2 + (j * dy)
                        Elec=simpsons1d(Electric,(-m+X)/2,(m - X)/2,Y,n)    #call function here to prevent having to call twice in intensity equation
                        intensity[i,j] = epsilon * c * (abs(Elec))**2   #equation to determine the intesity
                plt.imshow(intensity)    #plots intesity graph
                plt.colorbar()           #creates a scale next to the graph
                plt.show()
            
        
        elif MyInput5 == 'e':
            
            print("You chose to go back to the main option screen.")
            
        else:
            print("This is not a valid choice, please try again.")
        
        
    elif MyInput != 'q':
        print('This is not a valid choice, please try again.')
             
print('You have chosen to finish - goodbye.')
