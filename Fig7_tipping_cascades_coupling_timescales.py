# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 16:42:06 2025

@author: Paul

Plotting script to create Figure 7 in the paper "The role of coupling and
timescales for interacting tipping elements"
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import scipy.io as sp
from matplotlib import rc

fontsize = 10
rc('font', **{'size' : fontsize})
rc('text', usetex=True)

# Load in data
# mat_contents = sp.loadmat("Downstream_tipping_boundaries_linear_v2.mat")
# epsilon = mat_contents['epsilon'][0]
# b_track = mat_contents['b_track'][0][0]
# b_prop = mat_contents['b_prop'][0][0]
# b_DBUB = mat_contents['b_DBUB'][0]
# b_UBDB = mat_contents['b_UBDB'][0]
# b_UBf = mat_contents['b_UBf'][0]

mat_contents = sp.loadmat("Downstream_tipping_boundaries_linear_v3.mat")
epsilon = mat_contents['epsilon'][0]
b_track = mat_contents['b_track'][0][0]
b_prop = mat_contents['b_prop'][0][0]
b_on = mat_contents['b_on'][0][0]
b_off = mat_contents['b_off'][0]
b_UaDB = mat_contents['b_UaDB'][0]


# Plotting
fig, ax = plt.subplots(1,1,figsize=(7.007874,5.256))
# plt.plot(b_UBDB,epsilon,'k:')
# plt.plot([np.nanmax(b_UBDB),np.nanmax(b_UBDB)],[1e-3,1e2],'k')
# plt.plot(b_DBUB,epsilon,'k')
# plt.plot(b_UBf,epsilon,'k--')
# plt.plot([b_track,b_track],[1e-3,1e2],'k')
# plt.plot([b_prop,b_prop],[1e-3,1e2],'k--')

# plt.fill_between([0,b_track],1e-3,1e2,edgecolor='none',alpha=0.3)
# plt.fill_betweenx(epsilon,b_track*np.ones(len(epsilon)),b_UBDB,edgecolor='none',alpha=0.3)
# plt.fill_betweenx(epsilon,b_UBDB,b_DBUB,edgecolor='none',alpha=0.3)
# plt.fill_between(b_DBUB,1e-3*np.ones(len(b_DBUB)),epsilon,edgecolor='none',alpha=0.3)

plt.plot(b_UaDB,epsilon,'k')
plt.plot([b_on,b_on],[1e-3,1e2],'k--')
plt.plot(b_off,epsilon,'k--')
plt.plot([b_track,b_track],[1e-3,1e2],'k')
plt.plot([b_prop,b_prop],[1e-3,1e2],'k')

plt.fill_between([0,b_track],1e-3,1e2,edgecolor='none',alpha=0.3)
plt.fill_between([b_track,b_prop],1e-3,1e2,edgecolor='none',alpha=0.3)
plt.fill_betweenx(epsilon,b_prop*np.ones(len(epsilon)),b_UaDB,edgecolor='none',alpha=0.3)
plt.fill_between(b_UaDB,1e-3*np.ones(len(b_UaDB)),epsilon,edgecolor='none',alpha=0.3)

plt.xlim(0.1,10)
plt.ylim(1e-3,1e2)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$b$')
plt.ylabel('$\epsilon$')
plt.tight_layout()

plt.text(0.16,0.65,'UB',ha='center',va='center', transform=ax.transAxes)
plt.gca().add_patch(Rectangle((0.12,0.6),0.08,0.02,edgecolor='tab:blue',facecolor='tab:blue', transform=ax.transAxes))

plt.annotate("", xy=(0.52, 0.01), xytext=(0.3, 0.01), arrowprops=dict(arrowstyle="->"))
plt.text(0.195,0.2,'DaUB',ha='center',va='center', transform=ax.transAxes)
plt.gca().add_patch(Rectangle((0.13,0.15),0.08,0.02,edgecolor='tab:blue',facecolor='tab:blue', transform=ax.transAxes))
plt.gca().add_patch(Rectangle((0.22,0.12),0.04,0.02,edgecolor='tab:orange',facecolor='tab:orange', transform=ax.transAxes))

plt.text(0.6,0.3,'DwUB',ha='center',va='center', transform=ax.transAxes)
plt.gca().add_patch(Rectangle((0.56,0.25),0.08,0.02,edgecolor='tab:blue',facecolor='tab:blue', transform=ax.transAxes))
plt.gca().add_patch(Rectangle((0.58,0.22),0.04,0.02,edgecolor='tab:orange',facecolor='tab:orange', transform=ax.transAxes))

plt.text(0.57,0.8,'DoUB',ha='center',va='center', transform=ax.transAxes)
plt.gca().add_patch(Rectangle((0.5,0.75),0.08,0.02,edgecolor='tab:blue',facecolor='tab:blue', transform=ax.transAxes))
plt.gca().add_patch(Rectangle((0.54,0.72),0.1,0.02,edgecolor='tab:orange',facecolor='tab:orange', transform=ax.transAxes))

plt.text(0.895,0.3,'UaDB',ha='center',va='center', transform=ax.transAxes)
plt.gca().add_patch(Rectangle((0.88,0.25),0.08,0.02,edgecolor='tab:blue',facecolor='tab:blue', transform=ax.transAxes))
plt.gca().add_patch(Rectangle((0.83,0.22),0.04,0.02,edgecolor='tab:orange',facecolor='tab:orange', transform=ax.transAxes))

plt.text(0.89,0.95,'UwDB',ha='center',va='center', transform=ax.transAxes)
plt.gca().add_patch(Rectangle((0.85,0.9),0.08,0.02,edgecolor='tab:blue',facecolor='tab:blue', transform=ax.transAxes))
plt.gca().add_patch(Rectangle((0.83,0.87),0.12,0.02,edgecolor='tab:orange',facecolor='tab:orange', transform=ax.transAxes))

plt.text(0.895,0.74,'UoDB',ha='center',va='center', transform=ax.transAxes)
plt.gca().add_patch(Rectangle((0.815,0.69),0.08,0.02,edgecolor='tab:blue',facecolor='tab:blue', transform=ax.transAxes))
plt.gca().add_patch(Rectangle((0.775,0.66),0.08,0.02,edgecolor='tab:orange',facecolor='tab:orange', transform=ax.transAxes))

fig.subplots_adjust(top=0.962,bottom=0.104,left=0.12,right=0.965,hspace=0.2,wspace=0.2)
plt.text(-0.1,0.8,r"$\textbf{Decelerating}$" "\n" r"$\textbf{cascades}$",ha='center',va='center', rotation=90, transform=ax.transAxes)
plt.text(-0.1,0.3,r"$\textbf{Accelerating}$" "\n" r"$\textbf{cascades}$",ha='center',va='center', rotation=90, transform=ax.transAxes)

fig.savefig('../Figures/Linear_coupling/Cascade_tipping_b_eps_tipping_regimes_lin_v6.pdf', format='pdf', dpi=800)