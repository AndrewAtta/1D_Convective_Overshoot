#An attempt to implament the Forward Euler Method

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg



def BW_Euler(x_min = 0, x_max = 10, t_min = 0, t_max = 5, i_max = 100, n_max = 100, alpha = 1): #arbitrarily set starting values
    
    #define grid

    x = np.linspace(x_min, x_max, i_max) #begin by creating a grid of points along the x-axis
    dx = x[1] - x[0] #define the step size delta_x
    t = np.linspace(t_min, t_max, n_max) #define gride for t-axis
    dt = t[1] - t[0] #define step size delta_t based on number of steps
    
    #values within coefficient matrix
    
    a = - (alpha*dt)/(dx**2) #a for the matrix
    b = (1 + (2*alpha*dt)/dx**2) #b for the matrix
    
    #initial solution to test (in this case we take guassian)
    
    n = 0
    
    sigma = 1 #this is arbitrrily chosen
    
    u = (1/(2*sigma)) * np.exp(((-(x-5)**2) / (2*sigma**2)))
    u_exact = (1/(np.sqrt(2*sigma**2 + 4*alpha*n))) * np.exp(((-(x-5)**2) / (2*sigma**2 + 4*alpha*n)))
    
    u0 = u
    #u_exact= u
    
    print(u)
    
    #define tri-diagnoal coefficient matrix in banded matrix form
    
    coeff = np.zeros((3,n_max))
    coeff[0,:][1:n_max] = a
    coeff[2,:][0:n_max-1] = a 
    coeff[1,:] = b


    #boundry conditions
    
    
    
    while n < n_max:
        
        #numerical solution
        #u[0]= u[0] - a*0.1 # not sure what the boundry conditions here will be
        #u[-1]=u[-1] - a*0.1
        u_sol = scipy.linalg.solve_banded((1,1),coeff,u)
        
       
        #analytical solution
        
        #u_exact = (1/(2*sigma)) * np.exp(((-(x-5)**2) / (2*sigma**2)))
        u_exact = (1/np.sqrt((2*sigma**2 + 4*alpha*n ))) * np.exp(((-(x-5)**2) / (2*sigma**2 + 4*alpha*n)))
       
        
        plt.clf()
        plt.ylim(-0.5,1)
        plt.xlim(-3,13)
        plt.xlabel ('x')
        plt.ylabel ('u')
        c1, = plt.plot(x,u0,label='initial conditions')
        c2, = plt.plot(x,u_sol,label='numerical solution')
        c3, = plt.plot(x,u_exact,ls='dashed',label='analytic solution')
        plt.legend(handles=[c1,c2,c3])
        plt.draw()
        plt.show(block = False)
        plt.pause(0.05)
        print(u_sol)
        u = u_sol
        n = n+dt
        


    plt.clf()
    plt.ylim(-0.5,1)
    plt.xlim(-3,1)
    plt.xlabel ('x')
    plt.ylabel ('u')
    c1, = plt.plot(x,u0,label='initial conditions')
    c2, = plt.plot(x,u_sol,label='numerical solution')
    c3, = plt.plot(x,u_exact,ls='dashed',label='analytic solution')
    plt.legend(handles=[c1,c2,c3])
    plt.draw()
    plt.show(block = True)

BW_Euler()
