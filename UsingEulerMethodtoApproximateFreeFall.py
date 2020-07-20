#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 15:07:50 2018

@author: saskiajackson
"""

import numpy as np
import matplotlib.pyplot as plt

def plot(y,ylabel):  #a function to plot a graph
    plt.plot(T, y, linewidth=1)  #t is plotted in all the graphs whereas the yaxis changes so is left as a variable
    plt.xlabel('Time (s)')
    plt.ylabel(ylabel)  #the y axis label changes with the y axis
    
    
def multplot(y,legend,ylabel):   #a function to plot a graph with multiple plots and a legend
    plt.plot(T,y, label=legend)   #both the y axis and the legend are variables
    plt.xlabel('Time(s)')
    plt.ylabel(ylabel)    #the y axis label changes with the y axis
    plt.legend()
    
def plotmenu(Y,V):   #a function to call the menu system for plotting either a height or velocity graph
    MyInput2 = input('Enter a choice "y" to plot the distance graph, "v" to plot the velocity graph and "e" to exit to the main menu screen: ')
    if MyInput2 == 'y':
        plot(Y, 'Height (m)')   #plots a height-time graph
        plt.show()
    elif MyInput2 == 'v':
        plot(V, 'Velocity (m/s)')    #plots a velocity-time graph
        plt.show()
    elif MyInput2 == 'e':
        print("You chose to go back to the main options screen.")   #to exit inner menu system
    else:
        print("This is not a valid choice, please try again.")    #if any other keys are pressed other than the ones defined above
    
def eulers(t,y,v):   #a function for euler method 
    while y >=0:
        
        y += dt*v      #the series equation for height
        V.append(v)    #adds the calculated values for height, velocity or time to their respective array
        Y.append(y)
        T.append(t)
        
        v -= dt*(g + (k*abs(v)*v)/m)   #the series equation for velocity 
        #after the append lines in order to ensure that v=0 is plotted
        t+=dt   #increases the time in equal increments 
        
def analytical(t,y,v):     #a function for the analytical method 
    while y >=2:
        
        v = -np.sqrt((m*g)/k)*np.tanh(np.sqrt((k*g)/m)*t)     #the series equations for the analytical method 
        y = y0 - m*np.log((np.cosh(np.sqrt((k*g)/m)*t))**2)/(2*k)
        Va.append(v)   #adds the calculated values for height, velocity or time to their respective array
        Ya.append(y)
        T.append(t)
        
          
        t+=dt   #increases the time in equal increments
        
def modified(t,y,v):
    while y >=0:
        vmid = v -(dt*(g + (k*abs(v)*v)/m))/2     #the midpoint series equation
        v -= dt*(g + (k*abs(vmid)*vmid)/m)
        y += dt*vmid
        
        V.append(v)    #adds the calculated values for height, velocity or time to their respective array
        Y.append(y)
        T.append(t)
            
        t+=dt    #increases the time in equal increments
        
def densityeulers(t,y,v):
    while y >=0:
        p = p0*np.exp(-y/h)     #the varying density equation
        k = (C*p*A)/2
        vmid = v -(dt*(g + (k*abs(v)*v)/m))/2
        v -= dt*(g + (k*abs(vmid)*vmid)/m)
        y += dt*vmid
            
        V.append(v)    #adds the calculated values for height, velocity or time to their respective array
        Y.append(y)
        T.append(t)
            
        t+=dt    #increases the time in equal increments
 
MyInput = '0'
while MyInput != 'q':
    MyInput = input(''''Enter a choice:  "a" for the Euler method for a height of 1m,
                 "b" for an alternate analytical model
                 "c" for the modified Euler method
                 "d" for the programme adapted for Felix Baumgartner's jump with variable density
                 "e" for a comparison of the various graphs for different start times and CA/m 
                 "f" for a comparison of the various graphs for the euler method and analytical method or
                 "q" to quit: ''')
    print('You entered the choice: ', MyInput)
    if MyInput == 'a':
    
        C = 1.15     #sets the values for the set conditions of the free fall
        p0 = 1.2
        A = 0.7
        k = (C*p0*A)/2
        m = 120
        g = 9.81
        MyInput4 = input("Enter a value for delta t, the change in time between the iterations: ")   #allows a user inputted value for delta t
        
        try:
            dt = float(MyInput4)    #tests whether the input for delta t is a number
        except ValueError:   #if an error occurs then the user is informed 
            print("The value you entered for delta t is not a number, please try again.")
        else:    #if no error occurs then the program continues as normal
            y0 = 1E3
            V = []    #opens empty arrays to be appended to
            Y = []
            T = []
        
            y = y0     #sets the initial conditions for the while loop
            t = 0
            v = 0
            
            eulers(t,y,v)    #calls the function to run euler method
        
            plotmenu(Y,V)     #calls the function to plot the graphs
            
            
    elif MyInput == 'b':
        C = 1.15       #sets the values for the set conditions of the free fall
        p0 = 1.2
        A = 0.7
        k = (C*p0*A)/2
        m = 120
        g = 9.81
        MyInput4 = input("Enter a value for delta t, the change in time between the iterations: ")    #allows a user inputted value for delta t
        
        try:
            dt = float(MyInput4)    #tests whether the input for delta t is a number
        except ValueError:      #if an error occurs then the user is informed
            print("The value you entered for delta t is not a number, please try again.")
        else:          #if no error occurs then the program continues as normal
            y0 = 1E3
            Va = []     #opens empty arrays to be appended to
            Ya = []
            T = []
        
            y = y0     #sets the initial conditions for the while loop
            t = 0
            v = 0
            
            analytical(t,y,v)     #calls the function to run the analytical method
        
            plotmenu(Ya,Va)    #calls the function to plot the graphs
            
    elif MyInput == 'c':
        
        C = 1.15      #sets the values for the set conditions of the free fall
        p0 = 1.2
        A = 0.7
        k = (C*p0*A)/2
        m = 120
        g = 9.81
        MyInput4 = input("Enter a value for delta t, the change in time between the iterations: ")    #allows a user inputted value for delta t
        
        try:
            dt = float(MyInput4)    #tests whether the input for delta t is a number
        except ValueError:     #if an error occurs then the user is informed
            print("The value you entered for delta t is not a number, please try again.")
        else:           #if no error occurs then the program continues as normal
            y0 = (1*(10**3))
            V = []     #opens empty arrays to be appended to
            Y = []
            T = []     
        
            y = y0     #sets the initial conditions for the while loop
            t = 0
            v = 0
        
            modified(t,y,v)     #calls the function to run the modified euler method
        
            plotmenu(Y,V)     #calls the function to plot the graphs
        
    elif MyInput == 'd':
        
        C = 1.15       #sets the values for the set conditions of the free fall
        p0 = 1.2
        A = 0.7
        m = 120
        g = 9.81
        h = 7.64*(10**3)
        MyInput4 = input("Enter a value for delta t, the change in time between the iterations: ")     #allows a user inputted value for delta t
        
        try:
            dt = float(MyInput4)     #tests whether the input for delta t is a number
        except ValueError:     #if an error occurs then the user is informed
            print("The value you entered for delta t is not a number, please try again.")
        else:         #if no error occurs then the program continues as normal
            y0 = 39045
            V = []      #opens empty arrays to be appended to
            Y = []
            T = []
         
            y = y0      #sets the initial conditions for the while loop
            t = 0
            v = 0
        
            densityeulers(t,y,v)       #calls the function to run the modified euler method with varying density 
            
            
            plotmenu(Y,V)       #calls the function to plot the graphs
            
            print('The maximum velocity reached was ', abs(min(V)), 'm/s')   #prints the maximum value of the velocity 
        
        
    elif MyInput == 'e':
        p0 = 1.2    #sets the values for the set conditions of the free fall
        g = 9.81
        dt = 0.05
        h = 7.64E3
        
        y0 = 50000
        
        MyInput5 = input('Enter a choice "h" to vary the start height, "c" to vary CA/m and "e" to exit to the main menu screen: ')
        if MyInput5 == 'h':
            
            C = 1.15    #sets the values for the set conditions of the free fall
            A = 0.7
            m = 120
            
            for i in range(0,y0,10000):    #varying the start height from 0 to 50km, increasing by 10km each time
                V = []       #opens empty arrays to be appended to
                Y = []
                T = []
                
                t = 0       #sets the initial conditions for the while loop
                v = 0
                y = i
            
                densityeulers(t,y,v)     #calls the function to run the modified euler method with varying density 
                
                print('The maximum velocity reached for a start height of ',i/1000, 'km is ',abs(min(V)))    #prints the maximum velocity values in order to compare them
            
                multplot(V, "Height = %d km"%(i/1000), 'Velocity (m/s)')      #calls the function to plot the graphs
            plt.show()    #has to be inline with the for loop to ensure it adds the plots onto the same graph
            
        elif MyInput5 == 'c':
            for i in range(1,5,1):
                
                V = []      #opens empty arrays to be appended to
                Y = []
                T = []
            
                t = 0       #sets the initial conditions for the while loop
                v = 0
                a = i/100
                y = 39045
                
                while y >=0:
                    p = p0*np.exp(-y/h)      #the varying density equation
                    k = p/2    #changed the equation for k in order to allow CA/m to vary
                    vmid = v -(dt*(g + (k*a*abs(v)*v)))/2     #euler equations but with variable a instead of CA/m
                    v -= dt*(g + (k*a*abs(vmid)*vmid))
                    y += dt*vmid
            
                    V.append(v)
                    Y.append(y)
                    T.append(t)
            
                    t+=dt
                    
                print('The maximum velocity reached for a value of CA/m of ',i/100, 'is ',abs(min(V)))      #prints the maximum velocity values in order to compare them 
                
                multplot(V, "$C_d$A/m = %s"%(str(round(a,2))),'Velocity (m/s)')     #calls the function to plot the graphs
                plt.legend()
            plt.show()
            
        elif MyInput5 == 'e':
            print("You chose to go back to the main options screen.")
        
        else:
            print("This is not a valid choice, please try again.")
        
    elif MyInput == 'f':    #not a section of the exercise, but added in to create comparison graphs
        
        C = 1.15      #sets the values for the set conditions of the free fall
        p0 = 1.2
        A = 0.7
        k = (C*p0*A)/2
        m = 120
        g = 9.81
        dt = 0.05
        
        y0 = 1E3
        V = []      #opens empty arrays to be appended to
        Y = []
        T = []
        
        y = y0      #sets the initial conditions for the while loop
        t = 0
        v = 0
        
        
    
        eulers(t,y,v)     #calls the function to run the euler method
            
        Va = []     #resets the arrays to zero in order to allow another graph to be plotted
        Ya = []
        T = []
        
        y = y0      #sets the initial conditions for the while loop
        t = 0
        v = 0
        
        analytical(t,y,v)    #calls the function to run the analytical method 
            
        Ya = np.array(Ya)  #sets the arrays in order for them to be mathematically manipulated
        Y = np.array(Y)
        V = np.array(V)
        Va = np.array(Va)
        Z = Y - Ya      #finds the difference between the arrays
        D = V - Va
                    
        plot(Z, 'Difference in height (m)')     #calls the function to plot the graphs
        plt.show()
        plt.clf()
        
        plot(D, ' Difference in velocity (m/s)')     #calls the function to plot the graphs
        plt.show()

    
    elif MyInput != 'q':
        print('This is not a valid choice.')
                
print('You have chosen to finish - goodbye.')