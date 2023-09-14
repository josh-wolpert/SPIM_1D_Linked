#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 18:03:36 2023

Plot results of linked 1D river network simulations

@author: joshwolpert
"""

# Import modules
import numpy as np
import matplotlib.pyplot as plt

#%% Load data
dt=25

# Area add scenarios
st0_z_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st0_elevationadd_10000000.0.csv', delimiter=',')
st1_z_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st1_elevationadd_10000000.0.csv', delimiter=',')
st2_z_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st2_elevationadd_10000000.0.csv', delimiter=',')
st3_z_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st3_elevationadd_10000000.0.csv', delimiter=',')
st4_z_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st4_elevationadd_10000000.0.csv', delimiter=',')

st0_z_105 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st0_elevationadd_10500000.0.csv', delimiter=',')
st1_z_105 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st1_elevationadd_10500000.0.csv', delimiter=',')
st2_z_105 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st2_elevationadd_10500000.0.csv', delimiter=',')
st3_z_105 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st3_elevationadd_10500000.0.csv', delimiter=',')
st4_z_105 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st4_elevationadd_10500000.0.csv', delimiter=',')

st0_z_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st0_elevationadd_11000000.0.csv', delimiter=',')
st1_z_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st1_elevationadd_11000000.0.csv', delimiter=',')
st2_z_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st2_elevationadd_11000000.0.csv', delimiter=',')
st3_z_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st3_elevationadd_11000000.0.csv', delimiter=',')
st4_z_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st4_elevationadd_11000000.0.csv', delimiter=',')

st0_z_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st0_elevationadd_12000000.0.csv', delimiter=',')
st1_z_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st1_elevationadd_12000000.0.csv', delimiter=',')
st2_z_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st2_elevationadd_12000000.0.csv', delimiter=',')
st3_z_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st3_elevationadd_12000000.0.csv', delimiter=',')
st4_z_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m1/st4_elevationadd_12000000.0.csv', delimiter=',')


st0_e_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st0_erosionadd_10000000.0.csv', delimiter=',')*(1000000/dt)
st1_e_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st1_erosionadd_10000000.0.csv', delimiter=',')*(1000000/dt)
st2_e_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st2_erosionadd_10000000.0.csv', delimiter=',')*(1000000/dt)
st3_e_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st3_erosionadd_10000000.0.csv', delimiter=',')*(1000000/dt)
st4_e_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st4_erosionadd_10000000.0.csv', delimiter=',')*(1000000/dt)

st0_e_105 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st0_erosionadd_10500000.0.csv', delimiter=',')*(1000000/dt)
st1_e_105 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st1_erosionadd_10500000.0.csv', delimiter=',')*(1000000/dt)
st2_e_105 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st2_erosionadd_10500000.0.csv', delimiter=',')*(1000000/dt)
st3_e_105 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st3_erosionadd_10500000.0.csv', delimiter=',')*(1000000/dt)
st4_e_105 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st4_erosionadd_10500000.0.csv', delimiter=',')*(1000000/dt)

st0_e_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st0_erosionadd_11000000.0.csv', delimiter=',')*(1000000/dt)
st1_e_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st1_erosionadd_11000000.0.csv', delimiter=',')*(1000000/dt)
st2_e_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st2_erosionadd_11000000.0.csv', delimiter=',')*(1000000/dt)
st3_e_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st3_erosionadd_11000000.0.csv', delimiter=',')*(1000000/dt)
st4_e_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st4_erosionadd_11000000.0.csv', delimiter=',')*(1000000/dt)

st0_e_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st0_erosionadd_12000000.0.csv', delimiter=',')*(1000000/dt)
st1_e_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st1_erosionadd_12000000.0.csv', delimiter=',')*(1000000/dt)
st2_e_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st2_erosionadd_12000000.0.csv', delimiter=',')*(1000000/dt)
st3_e_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st3_erosionadd_12000000.0.csv', delimiter=',')*(1000000/dt)
st4_e_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st4_erosionadd_12000000.0.csv', delimiter=',')*(1000000/dt)

# Area loss scenarios
st0_z_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st0_elevationloss_10000000.0.csv', delimiter=',')
st1_z_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st1_elevationloss_10000000.0.csv', delimiter=',')
st2_z_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st2_elevationloss_10000000.0.csv', delimiter=',')
st3_z_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st3_elevationloss_10000000.0.csv', delimiter=',')
st4_z_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st4_elevationloss_10000000.0.csv', delimiter=',')

st0_z_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st0_elevationloss_11000000.0.csv', delimiter=',')
st1_z_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st1_elevationloss_11000000.0.csv', delimiter=',')
st2_z_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st2_elevationloss_11000000.0.csv', delimiter=',')
st3_z_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st3_elevationloss_11000000.0.csv', delimiter=',')
st4_z_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st4_elevationloss_11000000.0.csv', delimiter=',')

st0_z_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st0_elevationloss_12000000.0.csv', delimiter=',')
st1_z_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st1_elevationloss_12000000.0.csv', delimiter=',')
st2_z_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st2_elevationloss_12000000.0.csv', delimiter=',')
st3_z_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st3_elevationloss_12000000.0.csv', delimiter=',')
st4_z_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st4_elevationloss_12000000.0.csv', delimiter=',')

st0_e_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st0_erosionloss_10000000.0.csv', delimiter=',')*(1000000/dt)
st1_e_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st1_erosionloss_10000000.0.csv', delimiter=',')*(1000000/dt)
st2_e_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st2_erosionloss_10000000.0.csv', delimiter=',')*(1000000/dt)
st3_e_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st3_erosionloss_10000000.0.csv', delimiter=',')*(1000000/dt)
st4_e_10 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st4_erosionloss_10000000.0.csv', delimiter=',')*(1000000/dt)

st0_e_105 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st0_erosionloss_10500000.0.csv', delimiter=',')*(1000000/dt)
st1_e_105 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st1_erosionloss_10500000.0.csv', delimiter=',')*(1000000/dt)
st2_e_105 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st2_erosionloss_10500000.0.csv', delimiter=',')*(1000000/dt)
st3_e_105 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st3_erosionloss_10500000.0.csv', delimiter=',')*(1000000/dt)
st4_e_105 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st4_erosionloss_10500000.0.csv', delimiter=',')*(1000000/dt)

st0_e_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st0_erosionloss_11000000.0.csv', delimiter=',')*(1000000/dt)
st1_e_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st1_erosionloss_11000000.0.csv', delimiter=',')*(1000000/dt)
st2_e_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st2_erosionloss_11000000.0.csv', delimiter=',')*(1000000/dt)
st3_e_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st3_erosionloss_11000000.0.csv', delimiter=',')*(1000000/dt)
st4_e_11 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st4_erosionloss_11000000.0.csv', delimiter=',')*(1000000/dt)

st0_e_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st0_erosionloss_12000000.0.csv', delimiter=',')*(1000000/dt)
st1_e_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st1_erosionloss_12000000.0.csv', delimiter=',')*(1000000/dt)
st2_e_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st2_erosionloss_12000000.0.csv', delimiter=',')*(1000000/dt)
st3_e_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st3_erosionloss_12000000.0.csv', delimiter=',')*(1000000/dt)
st4_e_12 = np.genfromtxt('/Volumes/PhD_1/PrecipitationPhaser_runs/SPIM_1D_Output/m15/st4_erosionloss_12000000.0.csv', delimiter=',')*(1000000/dt)

x0 = np.linspace(0,20000,201)
x1 = np.linspace(0,3000,31)+6000
x2 = np.linspace(0,3000,31)+9500
x3 = np.linspace(0,3000,31)+13000
x4 = np.linspace(0,3000,31)+16500

# Color Scheme
c1 = np.array([0,51,102])/255
c2 = np.array([0,128,255])/255
c3 = np.array([0,255,255])/255


#%% Plot Stream Profiles for area loss scenario
fig1,ax1 = plt.subplots(figsize = (6,4))
stream0 = ax1.plot(x0[:-1],st0_z_10[:-1],linewidth=4,color=c1)
stream1 = ax1.plot(x1[:-1],st1_z_10[:-1],linewidth=4,color=c1)
stream2 = ax1.plot(x2[:-1],st2_z_10[:-1],linewidth=4,color=c1)
stream3 = ax1.plot(x3[:-1],st3_z_10[:-1],linewidth=4,color=c1)
stream4 = ax1.plot(x4[:-1],st4_z_10[:-1],linewidth=4,color=c1)

stream5 = ax1.plot(x0[:-1],st0_z_105[:-1],linewidth=4,color=c2)
stream6 = ax1.plot(x1[:-1],st1_z_105[:-1],linewidth=4,color=c2)
stream7 = ax1.plot(x2[:-1],st2_z_105[:-1],linewidth=4,color=c2)
stream8 = ax1.plot(x3[:-1],st3_z_105[:-1],linewidth=4,color=c2)
stream9 = ax1.plot(x4[:-1],st4_z_105[:-1],linewidth=4,color=c2)

stream10 = ax1.plot(x0[:-1],st0_z_12[:-1],linewidth=4,color=c3)
stream11 = ax1.plot(x1[:-1],st1_z_12[:-1],linewidth=4,color=c3)
stream12 = ax1.plot(x2[:-1],st2_z_12[:-1],linewidth=4,color=c3)
stream13 = ax1.plot(x3[:-1],st3_z_12[:-1],linewidth=4,color=c3)
stream14 = ax1.plot(x4[:-1],st4_z_12[:-1],linewidth=4,color=c3)

ax1.set_ylim(0,1000)
ax1.set_xlim([0,20000])
ax1.grid(False)
labels = [item.get_text() for item in ax1.get_xticklabels()]
labels = ax1.get_xticks().tolist()
labels = [round(item/1000) for item in labels]
ax1.set_xticklabels(labels)
ax1.set(xlabel='Distance (m)', ylabel='Elevation (km)')
ax1.set_facecolor((1, 1, 1))
ax1.spines['bottom'].set_color('0.5')
ax1.spines['top'].set_color('0.5')
ax1.spines['right'].set_color('0.5')
ax1.spines['left'].set_color('0.5')
#fig1.savefig('SPIM_1D_Output/profiles_arealoss.png', facecolor='white', edgecolor='none',dpi=500)


#%% Plot Erosion Rates for area loss scenario
fig1,ax1 = plt.subplots(figsize = (6,4))
stream0 = ax1.plot(x0[1:-1],st0_e_10[1:-1],linewidth=4,color=c1)
stream1 = ax1.plot(x1[1:-1],st1_e_10[1:-1],linewidth=4,color=c1)
stream2 = ax1.plot(x2[1:-1],st2_e_10[1:-1],linewidth=4,color=c1)
stream3 = ax1.plot(x3[1:-1],st3_e_10[1:-1],linewidth=4,color=c1)
stream4 = ax1.plot(x4[1:-1],st4_e_10[1:-1],linewidth=4,color=c1)

stream5 = ax1.plot(x0[1:-1],st0_e_105[1:-1],linewidth=4,color=c2)
stream6 = ax1.plot(x1[1:-1],st1_e_105[1:-1],linewidth=4,color=c2)
stream7 = ax1.plot(x2[1:-1],st2_e_105[1:-1],linewidth=4,color=c2)
stream8 = ax1.plot(x3[1:-1],st3_e_105[1:-1],linewidth=4,color=c2)
stream9 = ax1.plot(x4[1:-1],st4_e_105[1:-1],linewidth=4,color=c2)

stream10 = ax1.plot(x0[1:-1],st0_e_11[1:-1],linewidth=4,color=c3)
stream11 = ax1.plot(x1[1:-1],st1_e_11[1:-1],linewidth=4,color=c3)
stream12 = ax1.plot(x2[1:-1],st2_e_11[1:-1],linewidth=4,color=c3)
stream13 = ax1.plot(x3[1:-1],st3_e_11[1:-1],linewidth=4,color=c3)
stream14 = ax1.plot(x4[1:-1],st4_e_11[1:-1],linewidth=4,color=c3)

#ax1.set_ylim(160,205)
ax1.set_xlim([0,20000])
ax1.grid(False)
labels = [item.get_text() for item in ax1.get_xticklabels()]
labels = ax1.get_xticks().tolist()
labels = [round(item/1000,2) for item in labels]
ax1.set_xticklabels(labels)
ax1.set(xlabel='Distance (km)', ylabel='Erosion Rate (m My$r^{-1}$)')
ax1.set_facecolor((1, 1, 1))
ax1.spines['bottom'].set_color('0.5')
ax1.spines['top'].set_color('0.5')
ax1.spines['right'].set_color('0.5')
ax1.spines['left'].set_color('0.5')
fig1.savefig('SPIM_1D_Output/m15/erates_areaadd.png', facecolor='white', edgecolor='none',dpi=500)

