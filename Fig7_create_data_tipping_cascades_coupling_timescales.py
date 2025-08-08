# -*- coding: utf-8 -*-
"""
Created on Tue Aug  5 11:45:44 2025

@author: Paul

Script for simulation and plotting of deterministic trajectories for coupled
upstream and downstream systems. Coupling is through the localised state (overshoot)
"""

import numpy as np
from scipy.io import savemat

def f(x, forcing):
    """
    Generic double fold bifurcation model used for both upstream and downstream
    systems.

    Parameters
    ----------
    x : float64
        State variable for upstream or downstream system.
    forcing : float64
        External forcing to drive the system.

    Returns
    -------
    f : float64
        RHS of ODE.

    """  
    f = 3*x - x**3 + forcing    
    return (f)

def forcing_ramp(t, lamb_min, lamb_max, r):
    """
    Tanh ramp forcing for upstream system.

    Parameters
    ----------
    t : float64
        Time
    lamb_min : float64
        Initial level of forcing.
    lamb_max : float64
        Final level of forcing.
    r : float64
        Rate parameter of ramp.

    Returns
    -------
    lamb : float64
        Upstream ramp forcing profile.

    """
    lamb = lamb_min + (lamb_max-lamb_min)*(np.tanh(r*t)+1)/2   
    return (lamb)

def forcing_lin(x, xeq_start, a, b):
    """
    Linear coupling from upstream to downstream.

    Parameters
    ----------
    x : float64
        State variable of upstream system.
    xeq_start : float64
        Starting equiibrium of upstream system
    a : float64
        Initial level of forcing.
    b : float64
        Rate of linear increase.

    Returns
    -------
    lin : float64
        Linear coupling for downstream forcing profile.

    """   
    lin = a + b*(x-xeq_start)
    
    return (lin)


# Time parameters
tstart = -55                        # Start time
tend = 50                           # End time
dt = 0.0005                          # Spacing between time intervals
n = int((tend - tstart)/dt)         # Number of time intervals
t = np.linspace(tstart, tend, n+1)  # Time values

# Neighbourhood of alternative state threshold
w = 1.8

# Upstream forcing parameters
lamb_min = 0                # Start level of upstream forcing
lamb_max = 4                # End level of upstream forcing 
lamb_u = 2                  # Location of upper fold bifurcation
r = 0.05                    # Rate parameter of upstream forcing
xeq_start = -np.sqrt(3)     # Starting equilibrium for upstream system

# Downstream forcing parameters
a = 0               # Start level of downstream forcing
c = 2               # Inv. prop. to overshoot period of downstream forcing
d = 0.5             # Offset for coupling of downstream forcing during tipping of upstream 

# Initialise upstream variable
x_det_low = np.zeros(n+1)
x_det_low[0] = xeq_start 

for i in range(n):
    x_det_low[i+1] = x_det_low[i] + dt*f(x_det_low[i], forcing_ramp(t[i], lamb_min, lamb_max, r))      

# Identify time index for upstream onset
idx_low = np.argmin(np.abs(forcing_ramp(t, lamb_min, lamb_max, r)-lamb_u))

# Calculate tracking, simultaneous upstream offset & downstream onset, and simultaneous onset boundaries 
b_track = (lamb_u - a)/(x_det_low[-1]-xeq_start)
b_prop = (lamb_u - a)/(w-xeq_start)
b_on = (lamb_u - a)/(x_det_low[idx_low]-xeq_start) 

# Decreased time step and create high resolution upstream data
dt = 0.0001                          # Spacing between time intervals
n = int((tend - tstart)/dt)         # Number of time intervals
t = np.linspace(tstart, tend, n+1)  # Time values

idx_hi = np.argmin(np.abs(forcing_ramp(t, lamb_min, lamb_max, r)-lamb_u))

x_det_hi = np.zeros(n+1)
x_det_hi[0] = xeq_start 

for i in range(n):
    x_det_hi[i+1] = x_det_hi[i] + dt*f(x_det_hi[i], forcing_ramp(t[i], lamb_min, lamb_max, r)) 


# Set tolerance level for bisection method
err = 0.00001

# Array of different timescale separation parameters to consider
epsilon = np.concatenate((np.geomspace(0.001,10,24),np.geomspace(10,100,20)[1:]))

## Bisection method to determine coupling strength for simultaneous upstream onset & downstream offset
b_UaDB = np.zeros(len(epsilon))

# Loop over all timescale separation parameter values
for j in range(len(epsilon)):
    
    # If epsilon is small choose smaller time step and use corresponding upstream data
    if epsilon[j]<=4e-3:
        dt = 0.0001                          # Spacing between time intervals
        n = int((tend - tstart)/dt)         # Number of time intervals
        t = np.linspace(tstart, tend, n+1)  # Time values
        x_det = x_det_hi.copy()
        idx = idx_hi.copy()
    else:
        dt = 0.0005                          # Spacing between time intervals
        n = int((tend - tstart)/dt)         # Number of time intervals
        t = np.linspace(tstart, tend, n+1)  # Time values
        x_det = x_det_low.copy()
        idx = idx_low.copy()

    # Initial lower and upper bounds for coupling strength
    ba = b_on
    bb = 15
    
    # Initialise error level
    eps = np.abs(bb - ba)
    
    while eps>err:
        
        # Calculate mid point for coupling strength
        bc = (bb + ba)/2
        
        # Initialise variables     
        y_det_lin = np.zeros(n+1)
        speed = np.zeros(n+1)
    
        y_det_lin[0] = xeq_start        
        
        # Perform simulation until upstream onset or offstream offset (whichever occurs first)
        for i in range(idx):   
            y_det_lin[i+1] = y_det_lin[i] + dt*f(y_det_lin[i], forcing_lin(x_det[i], xeq_start, a, bc))/epsilon[j]
            
            if y_det_lin[i+1]>w:
                break
        
        # If downstream tipping finishes before upstream onset require weaker coupling strength
        if y_det_lin[i+1]>w:
            bb = bc
        else:
            ba = bc
        
        # Calculate new mid point
        bc = (bb + ba)/2
        
        # Calculate current difference between upper and lower bounds
        eps = np.abs(bb - ba)
    
    # Store coupling strength corresponding to simultaneous upstream onset & downstream offset
    b_UaDB[j] = bc


## Bisection method to determine coupling strength for simultaneous offset boundary
b_off = np.zeros(len(epsilon))

# Loop over all timescale separation parameter values
for j in range(len(epsilon)):
    
    # If epsilon is small choose smaller time step and use corresponding upstream data
    if epsilon[j]<=4e-3:
        dt = 0.0001                          # Spacing between time intervals
        n = int((tend - tstart)/dt)         # Number of time intervals
        t = np.linspace(tstart, tend, n+1)  # Time values
        x_det = x_det_hi.copy()
    else:
        dt = 0.0005                          # Spacing between time intervals
        n = int((tend - tstart)/dt)         # Number of time intervals
        t = np.linspace(tstart, tend, n+1)  # Time values
        x_det = x_det_low.copy()

    # Initial lower and upper bounds for coupling strength
    ba = b_track
    bb = 15
    
    # Initialise error level
    eps = np.abs(bb - ba)
    
    # While error is above tolerance perform bisection method
    while eps>err:
        
        # Calculate mid point for coupling strength
        bc = (bb + ba)/2
        
        # Initialise variables     
        y_det_lin = np.zeros(n+1)
    
        y_det_lin[0] = xeq_start        
        
        for i in range(n):   
            y_det_lin[i+1] = y_det_lin[i] + dt*f(y_det_lin[i], forcing_lin(x_det[i], xeq_start, a, bc))/epsilon[j]

            if (x_det[i+1]>w) or (y_det_lin[i+1]>w):
                break
            
        # If downstream tipping finishes first require weaker coupling strength
        if y_det_lin[i+1]>w:
            bb = bc
        else:
            ba = bc
        
        # Calculate new mid point
        bc = (bb + ba)/2
        
        # Calculate current difference between upper and lower bounds
        eps = np.abs(bb - ba)
    
    # Store coupling strength corresponding to simultaneous offset boundary
    b_off[j] = bc

    
# Save data
mdic = {"epsilon": epsilon, "b_track": b_track, "b_prop": b_prop, "b_on": b_on, "b_off": b_off, "b_UaDB": b_UaDB}

savemat("Downstream_tipping_boundaries_linear_v3.mat", mdic)