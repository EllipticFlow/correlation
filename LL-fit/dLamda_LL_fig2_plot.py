#!/usr/bin/env python
# coding: utf-8


import math
import torch
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfi
from torch.distributions.transforms import SoftmaxTransform
import torch.nn as nn

hbarc = 197.327

#define the |f(k)|^2=|1/(-1/a0+r0*k^2/2-i*k)|^2
def sfk(k,a,r):
    return 1/(k**2 + (-1/a + (k**2*r)/2)**2)

#define the Real part of f(k)
def Refk(k,a,r):
    return (k**2*r)/(2*(((k**2*r)/2 - 1/a)**2 + k**2)) - 1/(a*(((k**2*r)/2 - 1/a)**2 + k**2))

#define the Imaginary part of f(k)
def Imfk(k,a,r):
    return k/((1/a - (k**2*r)/2)**2 + k**2)

#define the F_0(r)
def F0(R,r):
    return 1-r/(2*math.sqrt(math.pi)*R)

#define the F_1(k)
def F1(k,R):
    return (math.sqrt(math.pi)*math.exp(-(2*k*R)**2)*erfi(2*k*R))/(2*(2*k*R))

#define the F_2(k)
def F2(k,R):
    return (1-math.exp(-(2*k*R)**2))/(2*k*R)

#define 1-state, eq(3) in paper arXiv:2005.05012
def LL1state(k,R,a,r):
    return sfk(k,a,r)*F0(R,r)/(2*R**2)+2*Refk(k,a,r)*F1(k,R)/(math.sqrt(math.pi)*R)-Imfk(k,a,r)*F2(k,R)/R

#define 1/3*D+2/3*Q+, eq(2) in paper arXiv:2005.05012
def CF(k,R,a1,r1,a2,r2):
    return 1+1/3*LL1state(k,R,a1,r1)+2/3*LL1state(k,R,a2,r2)

#define if only 1/3*D
def CF_D(k,R,a1,r1):
    return 1+1/3*LL1state(k,R,a1,r1)

#define if only 2/3*D
def CF_Q(k,R,a2,r2):
    return 1+2/3*LL1state(k,R,a2,r2)

#same as CF, just for consistency check
def CF_cob(k,R,a1,r1,a2,r2):
    return CF_D(k,R,a1,r1)+CF_Q(k,R,a2,r2)-1


# Define k values to plot CF
k_values = torch.linspace(0.2, 200, 100)


######################################################################
#####################################################################
# Compute CF values for each k value, fig2.L, 4SA
R = 2.5
a1 = 16.8 
r1 = 2.3 
a2 = -17.3
r2 = 3.6
CF_values_4SA = [CF(k/hbarc, R, a1, r1, a2, r2) for k in k_values]
CF_values_D_4SA = [CF_D(k/hbarc, R, a1, r1) for k in k_values]
CF_values_Q_4SA = [CF_Q(k/hbarc, R, a2, r2) for k in k_values]
CF_values_cob_4SA = [CF_cob(k/hbarc, R, a1, r1, a2, r2) for k in k_values]

# Compute CF values for each k value, fig2.L, 4Sf
R = 2.5
a1 = 16.8 
r1 = 2.3 
a2 = -10.8
r2 = 3.8
CF_values_4Sf = [CF(k/hbarc, R, a1, r1, a2, r2) for k in k_values]
CF_values_D_4Sf = [CF_D(k/hbarc, R, a1, r1) for k in k_values]
CF_values_Q_4Sf = [CF_Q(k/hbarc, R, a2, r2) for k in k_values]
CF_values_cob_4Sf = [CF_cob(k/hbarc, R, a1, r1, a2, r2) for k in k_values]

# Compute CF values for each k value, fig2.L, 4SE
R = 2.5
a1 = 16.8 
r1 = 2.3 
a2 = -7.6
r2 = 3.6
CF_values_4SE = [CF(k/hbarc, R, a1, r1, a2, r2) for k in k_values]
CF_values_D_4SE = [CF_D(k/hbarc, R, a1, r1) for k in k_values]
CF_values_Q_4SE = [CF_Q(k/hbarc, R, a2, r2) for k in k_values]
CF_values_cob_4SE = [CF_cob(k/hbarc, R, a1, r1, a2, r2) for k in k_values]


# Create a figure and a grid of subplots 
fig = plt.figure(figsize=(8, 6))
ax00 = plt.subplot2grid((3, 2), (0, 0))
ax10 = plt.subplot2grid((3, 2), (1, 0))
ax20 = plt.subplot2grid((3, 2), (2, 0))
ax01 = plt.subplot2grid((3, 2), (0, 1))
ax11 = plt.subplot2grid((3, 2), (1, 1))
ax21 = plt.subplot2grid((3, 2), (2, 1))

fig.suptitle("arXiv:2005.05012; Fig.2(L); R=2.5 fm")
#reduce the space between plot and title
fig.subplots_adjust(top=0.92) 

# Plot the data on each subplot

# Plot CF values as a function of k ###4SE -RED
ax00.plot(k_values, CF_values_4SE, color='#ff0000', label='2S+4SE')
ax00.plot(k_values, CF_values_D_4SE, color='#ff8800', label='1/3*2S')
ax00.plot(k_values, CF_values_Q_4SE, color='#a52a2a', label='2/3*4SE')
ax00.plot(k_values, CF_values_cob_4SE, color='black', linestyle=':',label='1/3+2/3')
ax00.legend(frameon=False,fontsize=8)
ax00.axhline(y=1, color='gray', linestyle='--') # Add dashed line at y=1
ax00.set_xlim(-5,155)
ax00.set_ylim(0,20)
# Plot CF values as a function of k
ax01.plot(k_values, CF_values_4SE, color='#ff0000', label='2S+4SE')
ax01.plot(k_values, CF_values_D_4SE, color='#ff8800', label='1/3*2S')
ax01.plot(k_values, CF_values_Q_4SE, color='#a52a2a', label='2/3*4SE')
ax01.plot(k_values, CF_values_cob_4SE, color='black', linestyle=':',label='1/3+2/3')
ax01.legend(frameon=False,fontsize=8)
ax01.text(5,1.1,"zoom")
ax01.axhline(y=1, color='gray', linestyle='--') # Add dashed line at y=1
ax01.set_xlim(-5,155)
ax01.set_ylim(0.9,1.2)


# Plot CF values as a function of k ###4SF -GREEN
ax10.plot(k_values, CF_values_4Sf, color='#008800', label='2S+4Sf')
ax10.plot(k_values, CF_values_D_4Sf, color='#00ff00', label='1/3*2S')
ax10.plot(k_values, CF_values_Q_4Sf, color='#338888', label='2/3*4Sf')
ax10.plot(k_values, CF_values_cob_4Sf, color='#ffe135', linestyle=':',label='1/3+2/3')
ax10.legend(frameon=False,fontsize=8)
ax10.axhline(y=1, color='gray', linestyle='--') # Add dashed line at y=1
ax10.set_xlim(-5,155)
ax10.set_ylim(0,20)
# Plot CF values as a function of k
ax11.plot(k_values, CF_values_4Sf, color='#008800', label='2S+4Sf')
ax11.plot(k_values, CF_values_D_4Sf, color='#00ff00', label='1/3*2S')
ax11.plot(k_values, CF_values_Q_4Sf, color='#338888', label='2/3*4Sf')
ax11.plot(k_values, CF_values_cob_4Sf, color='#ffe135', linestyle=':',label='1/3+2/3')
ax11.legend(frameon=False,fontsize=8)
ax11.text(5,1.1,"zoom")
ax11.axhline(y=1, color='gray', linestyle='--') # Add dashed line at y=1
ax11.set_xlim(-5,155)
ax11.set_ylim(0.9,1.2)

# Plot CF values as a function of k, ###4SA, Blue
ax20.plot(k_values, CF_values_4SA, color='#0066cc', label='2S+4SA')
ax20.plot(k_values, CF_values_D_4SA, color='#89cff0', label='1/3*2S')
ax20.plot(k_values, CF_values_Q_4SA, color='#800080', label='2/3*4SA')
ax20.plot(k_values, CF_values_cob_4SA, color='orange', linestyle=':',label='1/3+2/3')
ax20.legend(frameon=False,fontsize=8)
ax20.axhline(y=1, color='gray', linestyle='--') # Add dashed line at y=1
ax20.set_xlim(-5,155)
ax20.set_ylim(0,20)
# Plot CF values as a function of k,###4SA, Blue
ax21.plot(k_values, CF_values_4SA, color='#0066cc', label='2S+4SA')
ax21.plot(k_values, CF_values_D_4SA, color='#89cff0', label='1/3*2S')
ax21.plot(k_values, CF_values_Q_4SA, color='#800080', label='2/3*4SA')
ax21.plot(k_values, CF_values_cob_4SA, color='orange', linestyle=':',label='1/3+2/3')
ax21.legend(frameon=False,fontsize=8)
ax21.text(5,1.1,"zoom")
ax21.axhline(y=1, color='gray', linestyle='--') # Add dashed line at y=1
ax21.set_xlim(-5,155)
ax21.set_ylim(0.9,1.2)

# Save the plot

fig.text(0.5,0.04,"k (MeV/c)",ha='center')
fig.text(0.04,0.5,"C(k)",va='center', rotation='vertical')
fig.savefig('correlation_function_fig2_left.jpg', dpi=300, bbox_inches='tight')
#fig.show()


########################################################################
########################################################################
### Now plot the right figure of the figure 2
R = 2.5
a1 = 16.3 
r1 = 3.2 
a2 = -17.3
r2 = 3.6
CF_values_r4SA = [CF(k/hbarc, R, a1, r1, a2, r2) for k in k_values]
CF_values_D_r4SA = [CF_D(k/hbarc, R, a1, r1) for k in k_values]
CF_values_Q_r4SA = [CF_Q(k/hbarc, R, a2, r2) for k in k_values]
CF_values_cob_r4SA = [CF_cob(k/hbarc, R, a1, r1, a2, r2) for k in k_values]

# Compute CF values for each k value, fig2.L, 4Sf
R = 2.5
a1 = 16.3 
r1 = 3.2 
a2 = -10.8
r2 = 3.8
CF_values_r4Sf = [CF(k/hbarc, R, a1, r1, a2, r2) for k in k_values]
CF_values_D_r4Sf = [CF_D(k/hbarc, R, a1, r1) for k in k_values]
CF_values_Q_r4Sf = [CF_Q(k/hbarc, R, a2, r2) for k in k_values]
CF_values_cob_r4Sf = [CF_cob(k/hbarc, R, a1, r1, a2, r2) for k in k_values]

# Compute CF values for each k value, fig2.L, 4SE
R = 2.5
a1 = 16.3 
r1 = 3.2 
a2 = -7.6
r2 = 3.6
CF_values_r4SE = [CF(k/hbarc, R, a1, r1, a2, r2) for k in k_values]
CF_values_D_r4SE = [CF_D(k/hbarc, R, a1, r1) for k in k_values]
CF_values_Q_r4SE = [CF_Q(k/hbarc, R, a2, r2) for k in k_values]
CF_values_cob_r4SE = [CF_cob(k/hbarc, R, a1, r1, a2, r2) for k in k_values]


# Create a figure and a grid of subplots 
fig2 = plt.figure(figsize=(8, 6))
bx00 = plt.subplot2grid((3, 2), (0, 0))
bx10 = plt.subplot2grid((3, 2), (1, 0))
bx20 = plt.subplot2grid((3, 2), (2, 0))
bx01 = plt.subplot2grid((3, 2), (0, 1))
bx11 = plt.subplot2grid((3, 2), (1, 1))
bx21 = plt.subplot2grid((3, 2), (2, 1))

fig2.suptitle("arXiv:2005.05012; Fig.2(R); R=2.5 fm")
#reduce the space between plot and title
fig2.subplots_adjust(top=0.92) 

# Plot the data on each subplot

# Plot CF values as a function of k ###4SE -RED
bx00.plot(k_values, CF_values_r4SE, color='#ff0000', label='2S+4SE')
bx00.plot(k_values, CF_values_D_r4SE, color='#ff8800', label='1/3*2S')
bx00.plot(k_values, CF_values_Q_r4SE, color='#a52a2a', label='2/3*4SE')
bx00.plot(k_values, CF_values_cob_r4SE, color='black', linestyle=':',label='1/3+2/3')
bx00.legend(frameon=False,fontsize=8)
bx00.axhline(y=1, color='gray', linestyle='--') # Add dashed line at y=1
bx00.set_xlim(-5,155)
bx00.set_ylim(0,20)
# Plot CF values as a function of k
bx01.plot(k_values, CF_values_r4SE, color='#ff0000', label='2S+4SE')
bx01.plot(k_values, CF_values_D_r4SE, color='#ff8800', label='1/3*2S')
bx01.plot(k_values, CF_values_Q_r4SE, color='#a52a2a', label='2/3*4SE')
bx01.plot(k_values, CF_values_cob_r4SE, color='black', linestyle=':',label='1/3+2/3')
bx01.legend(frameon=False,fontsize=8)
bx01.text(5,1.1,"zoom")
bx01.axhline(y=1, color='gray', linestyle='--') # Add dashed line at y=1
bx01.set_xlim(-5,155)
bx01.set_ylim(0.9,1.2)


# Plot CF values as a function of k ###4SF -GREEN
bx10.plot(k_values, CF_values_r4Sf, color='#008800', label='2S+4Sf')
bx10.plot(k_values, CF_values_D_r4Sf, color='#00ff00', label='1/3*2S')
bx10.plot(k_values, CF_values_Q_r4Sf, color='#338888', label='2/3*4Sf')
bx10.plot(k_values, CF_values_cob_r4Sf, color='#ffe135', linestyle=':',label='1/3+2/3')
bx10.legend(frameon=False,fontsize=8)
bx10.axhline(y=1, color='gray', linestyle='--') # Add dashed line at y=1
bx10.set_xlim(-5,155)
bx10.set_ylim(0,20)
# Plot CF values as a function of k
bx11.plot(k_values, CF_values_r4Sf, color='#008800', label='2S+4Sf')
bx11.plot(k_values, CF_values_D_r4Sf, color='#00ff00', label='1/3*2S')
bx11.plot(k_values, CF_values_Q_r4Sf, color='#338888', label='2/3*4Sf')
bx11.plot(k_values, CF_values_cob_r4Sf, color='#ffe135', linestyle=':',label='1/3+2/3')
bx11.legend(frameon=False,fontsize=8)
bx11.text(5,1.1,"zoom")
bx11.axhline(y=1, color='gray', linestyle='--') # Add dashed line at y=1
bx11.set_xlim(-5,155)
bx11.set_ylim(0.9,1.2)

# Plot CF values as a function of k, ###4SA, Blue
bx20.plot(k_values, CF_values_r4SA, color='#0066cc', label='2S+4SA')
bx20.plot(k_values, CF_values_D_r4SA, color='#89cff0', label='1/3*2S')
bx20.plot(k_values, CF_values_Q_r4SA, color='#800080', label='2/3*4SA')
bx20.plot(k_values, CF_values_cob_r4SA, color='orange', linestyle=':',label='1/3+2/3')
bx20.legend(frameon=False,fontsize=8)
bx20.axhline(y=1, color='gray', linestyle='--') # Add dashed line at y=1
bx20.set_xlim(-5,155)
bx20.set_ylim(0,20)
# Plot CF values as a function of k,###4SA, Blue
bx21.plot(k_values, CF_values_r4SA, color='#0066cc', label='2S+4SA')
bx21.plot(k_values, CF_values_D_r4SA, color='#89cff0', label='1/3*2S')
bx21.plot(k_values, CF_values_Q_r4SA, color='#800080', label='2/3*4SA')
bx21.plot(k_values, CF_values_cob_r4SA, color='orange', linestyle=':',label='1/3+2/3')
bx21.legend(frameon=False,fontsize=8)
bx21.text(5,1.1,"zoom")
bx21.axhline(y=1, color='gray', linestyle='--') # Add dashed line at y=1
bx21.set_xlim(-5,155)
bx21.set_ylim(0.9,1.2)


fig2.text(0.5,0.04,"k (MeV/c)",ha='center')
fig2.text(0.04,0.5,"C(k)",va='center', rotation='vertical')
fig2.savefig('correlation_function_fig2_right.jpg', dpi=300, bbox_inches='tight')

