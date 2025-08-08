# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 11:44:57 2025

@author: Paul

Script to create Figure 8 in the paper "The role of coupling and
timescales for interacting tipping elements"
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import matplotlib as mpl
import matplotlib.colors as colors

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs
        
    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


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

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

fontsize = 10
rc('font', **{'size' : fontsize})
rc('text', usetex=True)

# Time parameters
tstart = -55                        # Start time
tend = 50                           # End time
dt = 0.0005                         # Spacing between time intervals
n = int((tend - tstart)/dt)         # Number of time intervals
t = np.linspace(tstart, tend, n+1)  # Time values
burn = 5                            # Burn time length to remove any transient behaviour
burn_pts = int(burn/dt)             # Number of points in burn period

# System parameter
epsilon = 0.05      # Timescale separation

# Upstream forcing parameters
lamb_min = 0                # Start level of upstream forcing
lamb_max = 4                # End level of upstream forcing 
r = 0.05                    # Rate parameter of upstream forcing
xeq_start = -np.sqrt(3)     # Starting equilibrium for upstream systeml

# Downstream forcing parameters
a = 0               # Start level of downstream forcing
c = 2               # Inv. prop. to overshoot period of downstream forcing
d = 0.5             # Offset for coupling of downstream forcing during tipping of upstream 

# Array of different coupling amplitudes b
bs = [1.8, 2, 2.2115194, 2.2115195, 8, 27.4, 29.3, 40]

# Initialise figure
fig2 = plt.figure(figsize=(7.007874,6))
# ax0 = fig2.add_axes([0.06,0.74,0.23,0.19])
ax0 = fig2.add_axes([0.73,0.08,0.2,0.2])
ax0.set_xlim(1.5,2.5)
ax0.set_ylim(0,4)
ax0.plot([2,2],[0,4],'k--',lw=1)
ax0.plot([0,4],[2,2],'k--',lw=1)
ax0.set_xlabel('$\lambda$')
ax0.set_ylabel('$M_2(x(\lambda))$')

# Set axes bounds 
lam_low_loc = -2.5
lam_high_loc = 5
x_high_loc = 4
y_high_loc = 4

# Locations of arrows to plot
larrowvectors = [[2], [2], [2,2.6], [2,2.6], [1.45,2.1], [1.45,2.2], [1.45,2.2], [1.45,2.2]]
xarrowvectors = [[0.2], [0.2], [0.2,1.9], [0.2,1.9], [0.9,0.4], [0.9,-0.2], [0.9,-0.2], [0.9,0]]
yarrowvectors = [[-2], [-2], [-2,-1.3], [-2,1.4], [0.9,2], [0.9,3], [0.9,3], [0.9,3.5]]

# List of colours for plotting
colours = ['Blues_r','Oranges_r','Greens_r','Reds_r','Purples_r','pink','bone','RdPu_r']

# Loop over the different coupling amplitudes b
for j in range(len(bs)):
    # Initialise variables 
    x_det = np.zeros(n+1)         
    y_det_loc = np.zeros(n+1)
    speed = np.zeros(n+1)
    
    # Initialise start points for variables
    x_det[0] = xeq_start
    y_det_loc[0] = xeq_start
    
    # Calculate deterministic trajectories
    for i in range(n):
        x_det[i+1] = x_det[i] + dt*f(x_det[i], forcing_ramp(t[i], lamb_min, lamb_max, r))
        y_det_loc[i+1] = y_det_loc[i] + dt*f(y_det_loc[i], forcing_loc(x_det[i], a, bs[j], c, d))/epsilon
        speed[i+1] = np.sqrt((x_det[i+1]-x_det[i])**2 + (y_det_loc[i+1]-y_det_loc[i])**2)

    # Set bounds (in log space) for colour intervals of different system speeds
    bounds = [-7,-4,-2,-1]
    cmap = truncate_colormap(plt.get_cmap(colours[j]), minval=0.2, maxval=0.7, n=100)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    # Plot panel (a)
    ax0.scatter(forcing_ramp(t, lamb_min, lamb_max, r),forcing_loc(x_det, a, bs[j], c, d),c=np.log10(speed),cmap=cmap,norm=norm,s=0.02,zorder=10, rasterized=True)
        
    # Initialise 3d subplot
    ax3 = fig2.add_subplot(3,3,j+1, projection='3d')
    ax3.view_init(elev=12., azim=-23)
    ax3.set_xlim(lam_low_loc,lam_high_loc)
    ax3.tick_params(axis='both', pad=-3)
    ax3.set_xticks([0,4])
    ax3.set_zticks([-2,0,2,4])
    ax3.set_xlabel('$\lambda$',labelpad=-7)
    ax3.set_ylabel('x',labelpad=-5)
    ax3.set_zlabel('y',labelpad=-10)
    
    ax3.xaxis.pane.fill = False
    ax3.yaxis.pane.fill = False
    ax3.zaxis.pane.fill = False
    
    # Now set color to white (or whatever is "invisible")
    ax3.xaxis.pane.set_edgecolor('w')
    ax3.yaxis.pane.set_edgecolor('w')
    ax3.zaxis.pane.set_edgecolor('w')

    ## Equilibria for localsied coupling
    yeq = np.linspace(-2,4,100001)
    yeq2 = np.linspace(-2,4,100001)
    xeq = np.arccosh(bs[j]/((yeq**3)-(3*yeq)-a))/c + d
    xeq2 = -np.arccosh(bs[j]/((yeq2**3)-(3*yeq2)-a))/c + d
    lambdaeq = (xeq**3) - (3*xeq)
    lambdaeq2 = (xeq2**3) - (3*xeq2)    
    
    lambdaeq[lambdaeq<lam_low_loc], xeq[lambdaeq<lam_low_loc], yeq[lambdaeq<lam_low_loc] = np.nan, np.nan, np.nan
    lambdaeq[lambdaeq>lam_high_loc], xeq[lambdaeq>lam_high_loc], yeq[lambdaeq>lam_high_loc] = np.nan, np.nan, np.nan
    lambdaeq[xeq>x_high_loc], xeq[xeq>x_high_loc], yeq[xeq>x_high_loc] = np.nan, np.nan, np.nan
    lambdaeq[yeq>y_high_loc], xeq[yeq>y_high_loc], yeq[yeq>y_high_loc] = np.nan, np.nan, np.nan
    lambdaeq2[lambdaeq2<lam_low_loc], xeq2[lambdaeq2<lam_low_loc], yeq2[lambdaeq2<lam_low_loc] = np.nan, np.nan, np.nan
    lambdaeq2[lambdaeq2>lam_high_loc], xeq2[lambdaeq2>lam_high_loc], yeq2[lambdaeq2>lam_high_loc] = np.nan, np.nan, np.nan
    lambdaeq2[xeq2>x_high_loc], xeq2[xeq2>x_high_loc], yeq2[xeq2>x_high_loc] = np.nan, np.nan, np.nan
    lambdaeq2[yeq2>y_high_loc], xeq2[yeq2>y_high_loc], yeq2[yeq2>y_high_loc] = np.nan, np.nan, np.nan
    
    # Plotting equilibria
    ax3.plot(lambdaeq,xeq,yeq,c='k',ls='--',lw=1)
    ax3.plot(lambdaeq[(3-3*(xeq**2)<0)&(3-3*(yeq**2)<0)],xeq[(3-3*(xeq**2)<0)&(3-3*(yeq**2)<0)],yeq[(3-3*(xeq**2)<0)&(3-3*(yeq**2)<0)],c='k',lw=1)
    ax3.plot(lambdaeq2,xeq2,yeq2,c='k',ls='--',lw=1)
    ax3.plot(lambdaeq2[(3-3*(xeq2**2)<0)&(3-3*(yeq2**2)<0)],xeq2[(3-3*(xeq2**2)<0)&(3-3*(yeq2**2)<0)],yeq2[(3-3*(xeq2**2)<0)&(3-3*(yeq2**2)<0)],c='k',lw=1)

    # Plotting system trajectory
    ax3.scatter(forcing_ramp(t[t>=-50], lamb_min, lamb_max, r),x_det[t>=-50],y_det_loc[t>=-50],c=np.log10(speed[t>=-50]),cmap=cmap,norm=norm,s=0.3,zorder=20, rasterized=True)#lw=2)
    
    # Add arrow(s) to trajectory
    larrow = larrowvectors[j]
    xarrow = xarrowvectors[j]
    yarrow = yarrowvectors[j]
    for k in range(len(larrow)):
        idx = np.nanargmin(np.abs(forcing_ramp(t, lamb_min, lamb_max, r) - larrow[k])+np.abs(x_det - xarrow[k])+np.abs(y_det_loc - yarrow[k]))
        arr = Arrow3D([forcing_ramp(t[idx], lamb_min, lamb_max, r), forcing_ramp(t[idx+1], lamb_min, lamb_max, r)], [x_det[idx], x_det[idx+1]], 
                    [y_det_loc[idx], y_det_loc[idx+1]], mutation_scale=12, 
                    lw=0, arrowstyle="simple", color=cmap(norm(np.log10(speed[idx]))),zorder=10)
        ax3.add_artist(arr)
    
fig2.subplots_adjust(top=1.0,bottom=0.0,left=0.0,right=0.965,hspace=0.0,wspace=0.0)
fig2.text(0.02,0.95,'\\textbf{(a)}')
fig2.text(0.34,0.95,'\\textbf{(b)}')
fig2.text(0.66,0.95,'\\textbf{(c)}')
fig2.text(0.02,0.62,'\\textbf{(d)}')
fig2.text(0.34,0.62,'\\textbf{(e)}')
fig2.text(0.66,0.62,'\\textbf{(f)}')
fig2.text(0.02,0.29,'\\textbf{(g)}')
fig2.text(0.34,0.29,'\\textbf{(h)}')
fig2.text(0.66,0.29,'\\textbf{(i)}')

fig2.text(0.28,0.95,'$b$ = '+str(bs[0]),ha='right')
fig2.text(0.6,0.95,'$b$ = '+str(bs[1]),ha='right')
fig2.text(0.92,0.95,'$b$ = '+str(bs[2]),ha='right')
fig2.text(0.28,0.62,'$b$ = '+str(bs[3]),ha='right')
fig2.text(0.6,0.62,'$b$ = '+str(bs[4]),ha='right')
fig2.text(0.92,0.62,'$b$ = '+str(bs[5]),ha='right')
fig2.text(0.28,0.29,'$b$ = '+str(bs[6]),ha='right')
fig2.text(0.6,0.29,'$b$ = '+str(bs[7]),ha='right')

fig2.text(0.945,0.07,'(a),(b),\n(c),(d)')
ax0.annotate("", xy=(2.36, 0.4), xytext=(2.57, 0.3), arrowprops=dict(arrowstyle="->"))
fig2.text(0.962,0.135,'(e)')
ax0.annotate("", xy=(2.4, 0.9), xytext=(2.66, 1.28), arrowprops=dict(arrowstyle="->"))
fig2.text(0.962,0.18,'(f)')
ax0.annotate("", xy=(2.5, 2.4), xytext=(2.66, 2.13), arrowprops=dict(arrowstyle="->"))
fig2.text(0.962,0.225,'(g)')
ax0.annotate("", xy=(2.5, 2.7), xytext=(2.66, 3.05), arrowprops=dict(arrowstyle="->"))
fig2.text(0.962,0.27,'(h)')
ax0.annotate("", xy=(2.5, 3.7), xytext=(2.66, 3.97), arrowprops=dict(arrowstyle="->"))

ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)

fig2.savefig('../Figures/Localised_coupling/Cascade_tipping_loc_coupling_lambda_M_3dbif_arrows_speed_v3.pdf', format='pdf', dpi=800)