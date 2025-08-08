# -*- coding: utf-8 -*-
"""
Created on Mon Aug  4 17:10:41 2025

@author: Paul

Script calculate boundary data for Figure 10 in the paper "The role of coupling and
timescales for interacting tipping elements"
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


def forcing_loc(x, a, b, c, d):
    """
    Localised state coupling from upstream to downstream.

    Parameters
    ----------
    x : float64
        State variable of upstream system.
    a : float64
        Initial level of forcing.
    b : float64
        Amplitude of forcing.
    c : float64
        DESCRIPTION.
    d : float64
        Offset coupling parameter.

    Returns
    -------
    loc : float64
        Localised state downstream overshoot forcing profile.

    """   
    loc = a + b/np.cosh(c*(x - d))
    
    return (loc)

# Upstream forcing parameters
lamb_min = 0                # Start level of upstream forcing
lamb_max = 4                # End level of upstream forcing 
lamb_u = 2                  # Location of upper fold bifurcation
r = 0.05                    # Rate parameter of upstream forcing

# Downstream forcing parameters
a = 0               # Start level of downstream forcing
c = 2               # Inv. prop. to overshoot period of downstream forcing
d = 0.5             # Offset for coupling of downstream forcing during tipping of upstream 

# Time parameters
tstart = -50                        # Start time
tend = 50                           # End time
dt = 0.01                           # Spacing between time intervals
n = int((tend - tstart)/dt)         # Number of time intervals
t = np.linspace(tstart, tend, n+1)  # Time values

# Set threshold for entering neighbourhood of alternative state
w = 1.8

# Identify time index at which upstream threshold is crossed 
idx = np.argmin(np.abs(forcing_ramp(t, lamb_min, lamb_max, r)-lamb_u))

# Initialise variables 
x_det = np.zeros(n+1)

x_coeff = [-1, 0, 3, forcing_ramp(tstart, lamb_min, lamb_max, r)]
x_det[0] = np.min(np.roots(x_coeff))

for i in range(n):
    x_det[i+1] = x_det[i] + dt*f(x_det[i], forcing_ramp(t[i], lamb_min, lamb_max, r))  

# Calculate boundary for simultaneous onset for upstream & downstream tipping
b_on = (lamb_u - a)*np.cosh(c*(x_det[idx]-d)) 

# Set tolerance level for bisection method
err = 0.001

# Array of different timescale separation parameters to consider
epsilon = np.concatenate((np.geomspace(0.001,10,24),np.geomspace(10,100,20)[1:]))

## Bisection method to determine coupling strength for boundary of cascading tipping
b_track = np.zeros(len(epsilon))

# Loop over all timescale separation parameter values
for j in range(len(epsilon)):
    
    # If epsilon is small choose smaller time step and extend simulation length
    if epsilon[j]/10>0.1:
        dt = 0.01
        tend = 500*epsilon[j]
    else:
        dt = epsilon[j]/10
        tend = 50

    n = int((tend - tstart)/dt)            # Number of time intervals
    t = np.linspace(tstart, tend, n+1)     # Time values
    
    # Initial lower and upper bounds for coupling strength
    ba = 0
    bb = 35    
    
    # Mid point for coupling strength
    bc = (bb + ba)/2   
    
    # Initialise error level
    eps = np.abs(bb - ba)
    
    # While error is above tolerance perform bisection method
    while eps>err:
        
        # Initialise variables 
        x_det = np.zeros(n+1)     
        y_det_loc = np.zeros(n+1)        
        
        x_coeff = [-1, 0, 3, forcing_ramp(tstart, lamb_min, lamb_max, r)]
        x_det[0] = np.min(np.roots(x_coeff))
        
        y_coeff = [-1, 0, 3, forcing_loc(x_det[0], a, bc, c, d)]
        y_det_loc[0] = np.min(np.roots(y_coeff))


        for i in range(n):
            x_det[i+1] = x_det[i] + dt*f(x_det[i], forcing_ramp(t[i], lamb_min, lamb_max, r))
            y_det_loc[i+1] = y_det_loc[i] + dt*f(y_det_loc[i], forcing_loc(x_det[i], a, bc, c, d))/epsilon[j]
            
            if y_det_loc[i+1]>1:
                break
            
        # If downstream system has tipped (note w = 1 here), require weaker coupling strength
        if y_det_loc[i+1] > 1:
            bb = bc
        else:
            ba = bc
            
        # Calculate new mid point
        bc = (bb + ba)/2
    
        # Calculate current difference between upper and lower bounds
        eps = np.abs(bb - ba)

    
    # Store coupling strength corresponding to boundary of cascading tipping    
    b_track[j] = bc



## Bisection method to determine coupling strength for simultaneous upstream onset & downstream offset
b_UaDB = np.zeros(len(epsilon))

# Loop over all timescale separation parameter values
for j in range(len(epsilon)):
    
    # If epsilon is small choose smaller time step and extend simulation length
    if epsilon[j]/10>0.1:
        dt = 0.01
        tend = 500*epsilon[j]
    else:
        dt = epsilon[j]/10
        tend = 50

    n = int((tend - tstart)/dt)         # Number of time intervals
    t = np.linspace(tstart, tend, n+1)  # Time values    

    # Identify time index for upstream onset 
    idx = np.argmin(np.abs(forcing_ramp(t, lamb_min, lamb_max, r)-lamb_u))
    
    # Initial lower and upper bounds for coupling strength
    ba = b_on
    bb = 90
    
    # Initialise error level
    eps = np.abs(bb - ba)
    
    while eps>err:
        
        # Calculate mid point for coupling strength
        bc = (bb + ba)/2
        
        # Initialise variables 
        x_det = np.zeros(n+1)     
        y_det_loc = np.zeros(n+1)        
        
        x_coeff = [-1, 0, 3, forcing_ramp(tstart, lamb_min, lamb_max, r)]
        x_det[0] = np.min(np.roots(x_coeff))
        
        y_coeff = [-1, 0, 3, forcing_loc(x_det[0], a, bc, c, d)]
        y_det_loc[0] = np.min(np.roots(y_coeff))

        # Perform simulation until upstream onset or offstream offset (whichever occurs first)
        for i in range(idx):
            x_det[i+1] = x_det[i] + dt*f(x_det[i], forcing_ramp(t[i], lamb_min, lamb_max, r))
            y_det_loc[i+1] = y_det_loc[i] + dt*f(y_det_loc[i], forcing_loc(x_det[i], a, bc, c, d))/epsilon[j]
            
            if y_det_loc[i+1]>w:
                break
        
        # If downstream tipping finishes before upstream onset require weaker coupling strength
        if y_det_loc[i+1]>w:
            bb = bc
        else:
            ba = bc
  
        # Calculate new mid point
        bc = (bb + ba)/2
        
        # Calculate current difference between upper and lower bounds
        eps = np.abs(bb - ba)
    
    # Store coupling strength corresponding to simultaneous upstream onset & downstream offset
    b_UaDB[j] = bc
    
    # Once coupling reaches a certain strength stop calculating and set all remaining to initial upper limit
    if bc > 85:
        b_UaDB[j+1:] = 85
        break


## Bisection method to determine coupling strength for simultaneous offset boundary
b_off = np.zeros(len(epsilon))

# Loop over all timescale separation parameter values
for j in range(len(epsilon)):
    
    # If epsilon is small choose smaller time step and extend simulation length
    if epsilon[j]/10>0.1:
        dt = 0.01
        tend = 500*epsilon[j]

    else:
        dt = epsilon[j]/10
        tend = 50
    
    n = int((tend - tstart)/dt)         # Number of time intervals
    t = np.linspace(tstart, tend, n+1)  # Time values
    
    # Initial lower and upper bounds for coupling strength
    ba = 0
    bb = 90
    
    # Initialise error level
    eps = np.abs(bb - ba)
    
    # While error is above tolerance perform bisection method
    while eps>err:
        
        # Calculate mid point for coupling strength
        bc = (bb + ba)/2
        
        # Initialise variables 
        x_det = np.zeros(n+1)     
        y_det_loc = np.zeros(n+1)        
        
        x_coeff = [-1, 0, 3, forcing_ramp(tstart, lamb_min, lamb_max, r)]
        x_det[0] = np.min(np.roots(x_coeff))
        
        y_coeff = [-1, 0, 3, forcing_loc(x_det[0], a, bc, c, d)]
        y_det_loc[0] = np.min(np.roots(y_coeff))


        for i in range(n):
            x_det[i+1] = x_det[i] + dt*f(x_det[i], forcing_ramp(t[i], lamb_min, lamb_max, r))
            y_det_loc[i+1] = y_det_loc[i] + dt*f(y_det_loc[i], forcing_loc(x_det[i], a, bc, c, d))/epsilon[j]

            if (x_det[i+1]>w) or (y_det_loc[i+1]>w):
                break
        
        # If downstream tipping finishes first require weaker coupling strength
        if y_det_loc[i+1]>w:
            bb = bc
        else:
            ba = bc
  
        # Calculate new mid point
        bc = (bb + ba)/2
        
        # Calculate current difference between upper and lower bounds
        eps = np.abs(bb - ba)
    
    # Store coupling strength corresponding to simultaneous offset boundary
    b_off[j] = bc
    
    # Once coupling reaches a certain strength stop calculating and set all remaining to initial upper limit
    if bc > 85:
        b_off[j+1:] = 85
        break


## Bisection method to determine theoretical coupling strength for boundary of cascading tipping
b_theory = np.zeros(len(epsilon))

# Loop over all timescale separation parameter values
for j in range(len(epsilon)):
    
    # Initial lower and upper bounds for coupling strength
    ba = 0
    bb = 35    
    
    # Mid point for coupling strength
    bc = (bb + ba)/2   
    
    # Initialise error level
    eps = np.abs(bb - ba)
    
    # While error is above tolerance perform bisection method
    while eps>err:
        
        # Initialise variables 
        x_det = np.zeros(n+1)            
        
        x_coeff = [-1, 0, 3, forcing_ramp(tstart, lamb_min, lamb_max, r)]
        x_det[0] = np.min(np.roots(x_coeff))
        
        y_coeff = [-1, 0, 3, forcing_loc(x_det[0], a, bc, c, d)]
        y_det_init = np.min(np.roots(y_coeff))

        for i in range(n):
            x_det[i+1] = x_det[i] + dt*f(x_det[i], forcing_ramp(t[i], lamb_min, lamb_max, r))
            
        forcing = forcing_loc(x_det, a, bc, c, d)
        
        # Calculate peak and time over threshold for downstream forcing 
        fmax = np.max(forcing)        
        tover = np.sum(forcing>2)*dt
        
        # If theory not satisfied require weaker coupling strength
        if (fmax-2)*(tover**2) > 4*(epsilon[j]**2)/3:
            bb = bc
        else:
            ba = bc
         
        # Calculate new mid point
        bc = (bb + ba)/2
        
        # Calculate current difference between upper and lower bounds
        eps = np.abs(bb - ba)   

    # Store theoretical coupling strength corresponding to boundary of cascading tipping
    b_theory[j] = bc

# Save data
mdic = {"epsilon": epsilon, "b_track": b_track, "b_on": b_on, "b_off": b_off, "b_UaDB": b_UaDB, "b_theory": b_theory}

savemat("Downstream_tipping_boundaries_theory_v5.mat", mdic)