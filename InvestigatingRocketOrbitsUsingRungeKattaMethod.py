#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 22:26:38 2018

@author: saskiajackson
"""

import matplotlib.pyplot as plt
def reset():    #resets the arrays so they can be appended to 
    
    Xaxis = []
    Yaxis = []
    Taxis = []
    GEaxis = []
    KEaxis = []
    TEaxis = []
    
    
    return(Xaxis, Yaxis, Taxis, GEaxis, KEaxis, TEaxis)

def dvx(X, Y):   #function to discribe the x component velocity differential for part a
    return ((-G*M*X)/((X**2 + Y**2)**(3/2)))

def dvy(X, Y):   #function to discribe the y component velocity differential for part a
    return ((-G*M*Y)/((X**2 + Y**2)**(3/2)))

def dvx_b(X,Y):    #function to discribe the x component velocity differential for part b of the exercise
    return (-((M*G*X)/(((X**2)+(Y**2))**(3/2)))-(Mm*G*X/((X**2+((ym-Y)**2))**(3/2))))

def dvy_b(X,Y):    #function to discribe the y component velocity differential for part b of the exercise
    return (-((M*G*Y)/(((X**2)+(Y**2))**(3/2)))+(Mm*G*(ym-Y)/((X**2+((ym-Y)**2))**(3/2))))
    

def runge(x,y,vx,vy,h):    #function to calculate the variables for the runge-kutta method
    k1x = vx 
    k1y = vy
    k1vx = dvx(x,y)    #calls the velocity equation functions 
    k1vy = dvy(x,y)
    
    k2x = vx + (h*k1vx)/2    #step size (h) is introduced into the equations
    k2y = vy + (h*k1vy)/2
    k2vx = dvx((x + (h*k1x)/2),(y + (h*k1y)/2))    
    k2vy = dvy((x + (h*k1x)/2),(y + (h*k1y)/2))
    
    k3x = vx + (h*k2vx)/2
    k3y = vy + (h*k2vy)/2
    k3vx = dvx((x + (h*k2x)/2),(y + (h*k2y)/2))
    k3vy = dvy((x + (h*k2x)/2),(y + (h*k2y)/2))
    
    k4x = vx + h*k3vx
    k4y = vy + h*k3vy
    k4vx = dvx((x + h*k3x),(y + h*k3y))
    k4vy = dvy((x + h*k3x),(y + h*k3y))
    return (k1x,k1y,k1vx,k1vy,k2x,k2y,k2vx,k2vy,k3x,k3y,k3vx,k3vy,k4x,k4y,k4vx,k4vy)   #returns the k values

def runge_b(x,y,vx,vy,h):    #function to calculate the variables for the runge-kutta method for part b
    k1x = vx 
    k1y = vy
    k1vx = dvx_b(x,y)  #calls the velocity equation functions
    k1vy = dvy_b(x,y)
    
    k2x = vx + (h*k1vx)/2      #stepsize is introduced
    k2y = vy + (h*k1vy)/2
    k2vx = dvx_b((x + (h*k1x)/2),(y + (h*k1y)/2))
    k2vy = dvy_b((x + (h*k1x)/2),(y + (h*k1y)/2))
    
    k3x = vx + (h*k2vx)/2
    k3y = vy + (h*k2vy)/2
    k3vx = dvx_b((x + (h*k2x)/2),(y + (h*k2y)/2))
    k3vy = dvy_b((x + (h*k2x)/2),(y + (h*k2y)/2))
    
    k4x = vx + h*k3vx
    k4y = vy + h*k3vy
    k4vx = dvx_b((x + h*k3x),(y + h*k3y))
    k4vy = dvy_b((x + h*k3x),(y + h*k3y))
    return (k1x,k1y,k1vx,k1vy,k2x,k2y,k2vx,k2vy,k3x,k3y,k3vx,k3vy,k4x,k4y,k4vx,k4vy)    #returns the k values 

def plot(distance,T):    #function to plot graphs of distance and energy 
    (Xaxis, Yaxis, Taxis, GEaxis, KEaxis, TEaxis) = reset()    #ensures the arrays are empty before appending anything to them 
    t = 0    #starts the time at 0
    h = 2    #step size
    
    if MyInput == 'a':   #if part a of the code is chosen, these values are set
        x = distance     #use of distance variable allows the start distance to be determined later
        y = 0
        vx = 0
        vy = (G*M/x)**(1/2)     #the start velocity of the rocket is determined by this equation to ensure a perfectly circular orbit
            
    elif MyInput == 'b':    #if part b of the code is chosen, these values are set
        
        vy = int(MyInput2)   
        x = distance      #use of distance variable allows the start distance to be determined later
        y = 0
        vx = 0
       
    elif MyInput == 'c':     #if part c of the code is chosen, these values are set
        x = 0
        y = distance      #use of distance variable allows the start distance to be determined later
        vx = -10.4815E3
        vy = 0
    
    for t in range(0,T):     #allows the time to be varied between 0 and the total time (set as a variable to allow it to be determined at a later point)
        
        if MyInput == 'a' or MyInput == 'b':    #if part a or b of the code is chosen, the first runge-kutta function is used
            
            (k1x,k1y,k1vx,k1vy,k2x,k2y,k2vx,k2vy,k3x,k3y,k3vx,k3vy,k4x,k4y,k4vx,k4vy) = runge(x,y,vx,vy,h)
        
        elif MyInput == 'c':      #if part c of the code is chosen, the second runge-kutta function is used
    
            (k1x,k1y,k1vx,k1vy,k2x,k2y,k2vx,k2vy,k3x,k3y,k3vx,k3vy,k4x,k4y,k4vx,k4vy) = runge_b(x,y,vx,vy,h)
            
        x += (h/6)*(k1x + 2*k2x + 2*k3x + k4x)       #sums the values of x,y,vx,vy and t calculated, causes them to increase in increments
        y += (h/6)*(k1y + 2*k2y + 2*k3y + k4y)
        vx += (h/6)*(k1vx + 2*k2vx + 2*k3vx + k4vx)
        vy += (h/6)*(k1vy + 2*k2vy + 2*k3vy + k4vy) 
        t += h
        
        if (x**2 + y**2)**(1/2) < 6371000:     #break function in order to ensure that the rocket doesnt crash into the Earth or the Moon
            break
        
        Xaxis.append(x)    #appends the distance and velocity values to the arrays
        Yaxis.append(y)
        Taxis.append(t)
        
        GE = (-G*M*m)/((x**2 + y**2)**(1/2))     #equation to determine the gravitational potential energy
        KE = (1/2)*(m*(vx**2 + vy**2))      #equation to determine the kinetic energy
        TE = GE + KE     #total energy is the sum of the gravitational and kinetic energies
         
        GEaxis.append(GE)      #appends the energy values to the array
        KEaxis.append(KE)
        TEaxis.append(TE)
        
    if MyInput == 'a' or MyInput == 'b':
        plt.plot(Xaxis, Yaxis)      #plots x against y to show a position plot
        plt.xlabel('Distance in x direction (m)')     #labels the x axis
        plt.ylabel('Distance in y direction (m)')     #labels the y axis
        plt.scatter(0,0,s=2E2,color='green')       #creates a green circle on the plot at position (0,0) with radius 200, to indicate the earth
        plt.gca().set_aspect('equal', adjustable= 'box')    #sets the axis of the graphs to be equal to allow the actual motion to be seen
        plt.show()     
        
        plt.plot(Taxis, GEaxis, label = 'Gravitational Energy')     #plots the gravitational energy against time and labels it
        plt.plot(Taxis, KEaxis, label = 'Kinetic Energy')     #plots the kinetic energy against time and labels it 
        plt.plot(Taxis, TEaxis, label = 'Total Energy')      #plots the total energy against time and labels it
        plt.xlabel('Time (s)')     #labels the x axis
        plt.ylabel('Energy (J)')     #labels the y axis
        plt.legend()    #shows the legend that discribes which line indicates which energy
        plt.show()
        
    elif MyInput == 'c':
        plt.plot(Xaxis, Yaxis)     #plots x against y to show a position plot
        plt.xlabel('Distance in x direction (m)')     #labels the x axis
        plt.ylabel('Distance in y direction (m)')     #labels the y axis
        plt.scatter(0,0,s=1.5E1,color='green')     #creates a green circle on the plot at position (0,0) with radius 15, to indicate the earth
        plt.scatter(0,384400E3, s=0.3E1, color = 'red')     #creates a red circle on the plot at position (0,384400E3) with radius 3, to indicate the moon
        plt.show()
        
        plt.plot(Taxis, GEaxis, label = 'Gravitational Energy')     #plots the gravitational energy against time and labels it
        plt.plot(Taxis, KEaxis, label = 'Kinetic Energy')      #plots the kinetic energy against time and labels it
        plt.plot(Taxis, TEaxis, label = 'Total Energy')       #plots the total energy against time and labels it
        plt.xlabel('Time (s)')     #labels the x axis
        plt.ylabel('Energy (J)')     #labels the y axis
        plt.legend()    #shows the legend that discribes which line indicates which energy
        plt.show()
    
        for i,j in enumerate(KEaxis):     #looks at the values for kinetic energy, in order to determine the period of the motion
            if round(j/(10**11),5) == 1.64812:     #if the value of kinetic energy is equal to 1.64812E11 then print the corresponding value of time
                print('The time taken for one period around the Earth and the Moon is: ', Taxis[i])    
         
        print('The closest the rocket gets to the Moon is: ', max(Yaxis) - 384400E3)   #finds the closest distance of the rocket of the Moon
        
G = 6.67E-11     #gravitational constant 
M = 5.9722E24    #mass of the Earth
m = 30000        #mass of the rocket

MyInput = '0'
while MyInput != 'q':
    MyInput = input('''Enter a choice:  "a" for a circular orbit around the Earth
                 "b" for an elliptical orbit around the earth
                 "c" for the motion of a rocket around the Earth and the Moon or
                 "q" to quit: ''')    #asks which part of the exercise the user would like to do
    
    if MyInput == 'a':
        print('You have chosen part (a); a circular orbit around the Earth')
        
        MyInput3 = input('Enter an integer value for the start distance from the Earth (for best results enter a number between 7,000,000 and 20,000,000): ')
              #allows user to input a value for the start distance
              
        try:    #checks if the inputted value is an integer
            X = int(MyInput3)
        except ValueError:
            print('''           ----------------------------------------------------
   The value entered for start distance is not an integer, please try again.
           ----------------------------------------------------''')
        else:
            T = 30000
            
            if X - 6371000 > 0:     #checks that the start distance is not smaller than the radius of the Earth
                
                plot(X,T)    #calls the plot function
                
            else:
                print('''           ----------------------------------------------------
   The value entered for start distance is smaller than the radius of the Earth
   please try again.
           ----------------------------------------------------''')
   
    elif MyInput == 'b':
        print('You have chosen part (b); an elliptical orbit around the Earth')
    
        MyInput2 = input('Enter an integer value for the start velocity in m/s (for most clear result enter a value between 6,000 and 9,000): ')
        MyInput3 = input('Enter an integer value for the start distance from the Earth (for best results enter a number between 7,000,000 and 20,000,000): ')
               #allows user to input a value of for the start distance and the start velocity
        
        try:    #checks if the inputted value is an integer
            vy = int(MyInput2)   
        except ValueError:
            print('''           ----------------------------------------------------
   The value entered for start velocity is not an integer, please try again.
           ----------------------------------------------------''')
        else:
            try:     #checks if the inputted value is an integer 
                X = int(MyInput3)
            except ValueError:
                print('''           ----------------------------------------------------
   The value entered for start distance is not an integer, please try again.
           ----------------------------------------------------''')
            else:
                T = 65000
                
                if X - 6371000 > 0:     #checks that the start distance is not smaller than the radius of the Earth
                    
                    plot(X,T)    #calls the plot function
                    
                else:
                    print('''           ----------------------------------------------------
   The value entered for start distance is smaller than the radius of the Earth, please try again.
           ----------------------------------------------------''')
           
    elif MyInput == 'c':
        print('You have chosen part (c); the motion of a rocket around the Earth and the Moon')
        
        ym = 384400E3    #distance from the Earth to the Moon
        Mm = 7.35E22     #mass of the Moon
        X = -71000E2     #start distance
        T = 500000       #maximum time of orbit
        
        plot(X,T)     #calls the plot function
        
    elif MyInput != 'q':
        print('''           ----------------------------------------------------
   This is not a valid option, please try again.
           ----------------------------------------------------''')
print('You have chosen to finish, goodbye')
        
