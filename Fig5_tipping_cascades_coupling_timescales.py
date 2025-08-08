# -*- coding: utf-8 -*-
"""
Created on Tue Aug  5 15:25:23 2025

@author: Paul

Script to create Figure 5 in the paper "The role of coupling and
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

# Downstream forcing parameter
a = 0               # Start level of downstream forcing

# Array of different coupling strengths b
bs = [0.495,0.52,1.25,3.468,6]

# Initialise figure
fig2 = plt.figure(figsize=(7.007874,4.5))
# ax0 = fig2.add_axes([0.06,0.58,0.27,0.32])
ax0 = fig2.add_axes([0.73,0.1,0.25,0.3])
ax0.set_xlim(0,4)
ax0.set_ylim(0,4)
ax0.plot([2,2],[0,4],'k--',lw=1)
ax0.plot([0,4],[2,2],'k--',lw=1)
ax0.set_xlabel(r'$\lambda$')
ax0.set_ylabel(r'$M_1(x(\lambda))$')

# Set axes bounds
lam_low_lin = -2.5
lam_high_lin = 5
x_high_lin = 4
y_low_lin = -2
y_high_lin = 4

# Locations of arrows to plot
larrowvectors = [[1,2.2], [1.2,2.2,3.2], [2,2.3,2.2], [1.3,2,2.2], [1.1,1.45,2.2]]
xarrowvectors = [[-1.7,0.6], [-1.7,0.6,0.9], [-1.7,0.9,1.3], [-1.5,-1.2,0.9], [-1.4,0.9,0.7]]
yarrowvectors = [[-1.7,-1.8], [-1.7,-1.8,0.9], [-1.7,0.9,2], [-1.5,0.9,2], [-1.4,0.9,2]]

# List of colours for plotting
colours = ['Blues_r','Oranges_r','Greens_r','Reds_r','Purples_r','pink','bone','RdPu_r']

# Loop over the different coupling strengths b
for j in range(len(bs)):
    # Initialise variables 
    x_det = np.zeros(n+1)     
    y_det_lin = np.zeros(n+1)
    speed = np.zeros(n)
    
    # Initialise start points for variables
    x_det[0] = xeq_start
    y_det_lin[0] = xeq_start
    
    # Calculate deterministic trajectories
    for i in range(n):
        x_det[i+1] = x_det[i] + dt*f(x_det[i], forcing_ramp(t[i], lamb_min, lamb_max, r))
        y_det_lin[i+1] = y_det_lin[i] + dt*f(y_det_lin[i], forcing_lin(x_det[i], xeq_start, a, bs[j]))/epsilon
        speed[i] = np.sqrt((x_det[i+1]-x_det[i])**2 + (y_det_lin[i+1]-y_det_lin[i])**2)
    
    # Set bounds (in log space) for colour intervals of different system speeds
    bounds = [-7,-4,-2,-1]
    cmap = truncate_colormap(plt.get_cmap(colours[j]), minval=0.2, maxval=0.7, n=100)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    # Plot panel (a)
    ax0.scatter(forcing_ramp(t[1:], lamb_min, lamb_max, r),forcing_lin(x_det[1:], xeq_start, a, bs[j]),c=np.log10(speed),cmap=cmap,norm=norm,s=0.02,zorder=10, rasterized=True)

    # Initialise 3d subplot
    ax2 = fig2.add_subplot(2,3,j+1, projection='3d')
    plt.rcParams['grid.color'] = [0.1,0.1,0.1,0.2]
    # ax2.view_init(elev=-155., azim=155)
    ax2.view_init(elev=12., azim=-23)
    ax2.set_xlim(lam_low_lin,lam_high_lin)
    ax2.set_zlim(y_low_lin,y_high_lin)
    ax2.tick_params(axis='both', pad=-3)
    ax2.set_xlabel(r'$\lambda$',labelpad=-5)
    ax2.set_ylabel('x',labelpad=-5)
    ax2.set_zlabel('y',labelpad=-10)
    
    ax2.xaxis.pane.fill = False
    ax2.yaxis.pane.fill = False
    ax2.zaxis.pane.fill = False
    
    # Now set color to white (or whatever is "invisible")
    ax2.xaxis.pane.set_edgecolor('w')
    ax2.yaxis.pane.set_edgecolor('w')
    ax2.zaxis.pane.set_edgecolor('w')

    ## Equilibria for linear coupling
    yeq = np.linspace(-2,5,1001)
    xeq = ((yeq**3)-(3*yeq)-a)/bs[j] + xeq_start
    lambdaeq = (xeq**3) - (3*xeq)    
    
    lambdaeq[lambdaeq<lam_low_lin], xeq[lambdaeq<lam_low_lin], yeq[lambdaeq<lam_low_lin] = np.nan, np.nan, np.nan
    lambdaeq[lambdaeq>lam_high_lin], xeq[lambdaeq>lam_high_lin], yeq[lambdaeq>lam_high_lin] = np.nan, np.nan, np.nan
    lambdaeq[xeq>x_high_lin], xeq[xeq>x_high_lin], yeq[xeq>x_high_lin] = np.nan, np.nan, np.nan
    lambdaeq[yeq>y_high_lin], xeq[yeq>y_high_lin], yeq[yeq>y_high_lin] = np.nan, np.nan, np.nan
    
    # Plotting equilibria
    ax2.plot(lambdaeq,xeq,yeq,c='k',ls='--',lw=1)
    ax2.plot(lambdaeq[(3-3*(xeq**2)<0)&(3-3*(yeq**2)<0)&(xeq>0)&(yeq>0)],xeq[(3-3*(xeq**2)<0)&(3-3*(yeq**2)<0)&(xeq>0)&(yeq>0)],yeq[(3-3*(xeq**2)<0)&(3-3*(yeq**2)<0)&(xeq>0)&(yeq>0)],c='k',lw=1)
    ax2.plot(lambdaeq[(3-3*(xeq**2)<0)&(3-3*(yeq**2)<0)&(xeq>0)&(yeq<=0)],xeq[(3-3*(xeq**2)<0)&(3-3*(yeq**2)<0)&(xeq>0)&(yeq<=0)],yeq[(3-3*(xeq**2)<0)&(3-3*(yeq**2)<0)&(xeq>0)&(yeq<=0)],c='k',lw=1)
    ax2.plot(lambdaeq[(3-3*(xeq**2)<0)&(3-3*(yeq**2)<0)&(xeq<=0)],xeq[(3-3*(xeq**2)<0)&(3-3*(yeq**2)<0)&(xeq<=0)],yeq[(3-3*(xeq**2)<0)&(3-3*(yeq**2)<0)&(xeq<=0)],c='k',lw=1)

    # Plotting system trajectory
    ax2.scatter(forcing_ramp(t[1:], lamb_min, lamb_max, r),x_det[1:],y_det_lin[1:],c=np.log10(speed),cmap=cmap,norm=norm,s=0.3,zorder=20, rasterized=True)#lw=2)

    # Add arrow(s) to trajectory
    larrow = larrowvectors[j]
    xarrow = xarrowvectors[j]
    yarrow = yarrowvectors[j]
    for k in range(len(larrow)):
        idx = np.nanargmin(np.abs(forcing_ramp(t, lamb_min, lamb_max, r) - larrow[k])+np.abs(x_det - xarrow[k])+np.abs(y_det_lin - yarrow[k]))
        arr = Arrow3D([forcing_ramp(t[idx], lamb_min, lamb_max, r), forcing_ramp(t[idx+1], lamb_min, lamb_max, r)], [x_det[idx], x_det[idx+1]], 
                    [y_det_lin[idx], y_det_lin[idx+1]], mutation_scale=14, 
                    lw=0, arrowstyle="simple", color=cmap(norm(np.log10(speed[idx]))),zorder=10)
        ax2.add_artist(arr)
    
# fig2.subplots_adjust(top=1.0,bottom=0.02,left=0.05,right=0.99,hspace=0.0,wspace=0.2)
fig2.subplots_adjust(top=1.0,bottom=0.0,left=0.0,right=0.965,hspace=0.0,wspace=0.1)
fig2.text(0.05,0.93,'\\textbf{(a)}')
fig2.text(0.38,0.93,'\\textbf{(b)}')
fig2.text(0.71,0.93,'\\textbf{(c)}')
fig2.text(0.05,0.43,'\\textbf{(d)}')
fig2.text(0.38,0.43,'\\textbf{(e)}')
fig2.text(0.71,0.43,'\\textbf{(f)}')

fig2.text(0.22,0.93,'$b$ = '+str(bs[0]))
fig2.text(0.55,0.93,'$b$ = '+str(bs[1]))
fig2.text(0.88,0.93,'$b$ = '+str(bs[2]))
fig2.text(0.22,0.43,'$b$ = '+str(bs[3]))
fig2.text(0.55,0.43,'$b$ = '+str(bs[4]))

fig2.text(0.94,0.15,'(a)')
ax0.annotate("", xy=(3.54, 1.85), xytext=(3.54, 1.1), arrowprops=dict(arrowstyle="->"))
fig2.text(0.94,0.33,'(b)')
ax0.annotate("", xy=(3.54, 2.15), xytext=(3.54, 2.9), arrowprops=dict(arrowstyle="->"))
fig2.text(0.94,0.38,'(c)')
ax0.annotate("", xy=(2.45, 3.85), xytext=(3.25, 3.85), arrowprops=dict(arrowstyle="->"))
fig2.text(0.77,0.3,'(d)')
ax0.annotate("", xy=(2.2, 2.8), xytext=(1.1, 2.8), arrowprops=dict(arrowstyle="->"))
fig2.text(0.77,0.35,'(e)')
ax0.annotate("", xy=(1.9, 3.48), xytext=(1.1, 3.48), arrowprops=dict(arrowstyle="->"))

ax0.spines['right'].set_visible(False)
ax0.spines['top'].set_visible(False)

# fig2.savefig('../Figures/Linear_coupling/Cascade_tipping_lin_coupling_lambda_M_3dbif_arrows_speed_v2.pdf', format='pdf', dpi=800)
