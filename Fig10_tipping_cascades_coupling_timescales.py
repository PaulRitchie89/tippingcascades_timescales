# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 17:23:32 2025

@author: Paul

Plotting script to create Figure 10 in the paper "The role of coupling and
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
mat_contents = sp.loadmat("Downstream_tipping_boundaries_theory_v5.mat")
epsilon = mat_contents['epsilon'][0]
b_track = mat_contents['b_track'][0]
b_on = mat_contents['b_on'][0][0]
b_off = mat_contents['b_off'][0]
b_UaDB = mat_contents['b_UaDB'][0]
b_theory = mat_contents['b_theory'][0]

# Plotting
fig, ax = plt.subplots(1,1,figsize=(7.007874,5.6))

plt.plot(b_UaDB,epsilon,'k')
plt.plot([b_on,b_on],[1e-3,1e2],'k--')
plt.plot(b_off,epsilon,'k--')
plt.plot(b_track,epsilon,'k')
plt.plot([2,2],[1e-3,1e2],'k:',alpha=0.5)
plt.plot([b_track[-1],b_track[-1]],[1e-3,1e2],'k:',alpha=0.5)

plt.fill_betweenx(epsilon,0,b_track,edgecolor='none',alpha=0.3)
plt.fill_betweenx(epsilon,b_track,b_UaDB,edgecolor='none',alpha=0.3,color='tab:green')
plt.fill_between(b_UaDB,1e-3*np.ones(len(b_UaDB)),epsilon,edgecolor='none',alpha=0.3,color='tab:red')

plt.plot(b_theory[epsilon<1],epsilon[epsilon<1],'r--')

plt.text(2.8,0.05,'Theory',ha='center',rotation=60,c='r')

plt.xlabel('$b$')
plt.ylabel('$\epsilon$')

plt.xlim(1,80)
plt.ylim(1e-3,1e2)

plt.xscale('log')
plt.yscale('log')


ax.annotate("", xy=(np.log10(2)/np.log10(80),1.03), xytext=(np.log10(b_track[-1])/np.log10(80), 1.03), xycoords='axes fraction', arrowprops=dict(arrowstyle="<->"))
plt.text((np.log10(2) +(np.log10(b_track[-1])-np.log10(2))/2)/np.log10(80),1.05,'Overshoot',ha='center',c='k', transform=ax.transAxes)


plt.text(0.28,0.8,'UB',ha='center',va='center', transform=ax.transAxes)
plt.gca().add_patch(Rectangle((0.25,0.75),0.06,0.02,edgecolor='tab:blue',facecolor='tab:blue', transform=ax.transAxes))

plt.text(0.49,0.35,'DwUB',ha='center',va='center', transform=ax.transAxes)
plt.gca().add_patch(Rectangle((0.46,0.3),0.06,0.02,edgecolor='tab:blue',facecolor='tab:blue', transform=ax.transAxes))
plt.gca().add_patch(Rectangle((0.475,0.27),0.03,0.02,edgecolor='tab:orange',facecolor='tab:orange', transform=ax.transAxes))

plt.text(0.64,0.72,'DoUB',ha='center',va='center', transform=ax.transAxes)
plt.gca().add_patch(Rectangle((0.66,0.78),0.06,0.02,edgecolor='tab:blue',facecolor='tab:blue', transform=ax.transAxes))
plt.gca().add_patch(Rectangle((0.69,0.75),0.06,0.02,edgecolor='tab:orange',facecolor='tab:orange', transform=ax.transAxes))

plt.text(0.88,0.3,'UaDB',ha='center',va='center', transform=ax.transAxes)
plt.gca().add_patch(Rectangle((0.87,0.25),0.06,0.02,edgecolor='tab:blue',facecolor='tab:blue', transform=ax.transAxes))
plt.gca().add_patch(Rectangle((0.83,0.22),0.03,0.02,edgecolor='tab:orange',facecolor='tab:orange', transform=ax.transAxes))

plt.text(0.875,0.96,'UwDB',ha='center',va='center', transform=ax.transAxes)
plt.gca().add_patch(Rectangle((0.845,0.91),0.06,0.02,edgecolor='tab:blue',facecolor='tab:blue', transform=ax.transAxes))
plt.gca().add_patch(Rectangle((0.83,0.88),0.09,0.02,edgecolor='tab:orange',facecolor='tab:orange', transform=ax.transAxes))

plt.text(0.87,0.74,'UoDB',ha='center',va='center', transform=ax.transAxes)
plt.gca().add_patch(Rectangle((0.8,0.69),0.06,0.02,edgecolor='tab:blue',facecolor='tab:blue', transform=ax.transAxes))
plt.gca().add_patch(Rectangle((0.76,0.66),0.06,0.02,edgecolor='tab:orange',facecolor='tab:orange', transform=ax.transAxes))

fig.subplots_adjust(top=0.962,bottom=0.104,left=0.12,right=0.965,hspace=0.2,wspace=0.2)
plt.text(-0.1,0.8,r"$\textbf{Decelerating}$" "\n" r"$\textbf{cascades}$",ha='center',va='center', rotation=90, transform=ax.transAxes)
plt.text(-0.1,0.3,r"$\textbf{Accelerating}$" "\n" r"$\textbf{cascades}$",ha='center',va='center', rotation=90, transform=ax.transAxes)

plt.tight_layout()
plt.subplots_adjust(top=0.93)

# fig.savefig('../Figures/Localised_coupling/Cascade_tipping_b_eps_tipping_regimes_loc_v6.pdf', format='pdf', dpi=800)