# -*- coding: utf-8 -*-
"""
Plot contraints of JH75 model
"""

import os
cd=os.getcwd()
import sys
sys.path.append('C:/Users/SLETIZIA/OneDrive - NREL/Desktop/PostDoc/utils')
import xarray as xr
import numpy as np
import utils as utl
from matplotlib import pyplot as plt
import warnings
import matplotlib
import glob

warnings.filterwarnings('ignore')
plt.close('all')

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.size'] = 14

#%% Inputs
source=os.path.join(cd,'data/old/DS*nc')
L_lim=[100,10000]#[m] limits of L
H_L_max=0.05#max slope 
H_l_max=1#max height over inner layer thickness
l_pblh_max=0.2#max inner layer height over pblh
pblh=1000#[m] height of PBL

#domain limits
xlim=[-250,1500]#[m]
ylim=[-875,875]#[m]

#%% Functions
def wrap(x):
    return np.append(x,x[0])

#%% Initalization
files=glob.glob(source)

#zeroing
wd=[]

H_all=[]
L_all=[]
l_all=[]

L_med=[]
L_low=[]
L_top=[]

H_L_med=[]
H_L_low=[]
H_L_top=[]

H_l_med=[]
H_l_low=[]
H_l_top=[]

l_med=[]
l_low=[]
l_top=[]

#%% Main
ctr=0
for f in files:
    Data=xr.open_dataset(f)
    wd=np.append(wd,np.float64(f[-6:-3]))
    Data=Data.where(Data['x']>xlim[0],drop=True).where(Data['x']<xlim[1],drop=True).where(Data['y']>ylim[0],drop=True).where(Data['y']<ylim[1],drop=True)
    
    H_all=np.append(H_all,Data.H.values)
    L_all=np.append(L_all,Data.L.values)
    l_all=np.append(l_all,Data.l.values)
    
    L_med=np.append(L_med,np.nanpercentile(Data.L.values,50))
    L_top=np.append(L_top,np.nanpercentile(Data.L.values,75))
    L_low=np.append(L_low,np.nanpercentile(Data.L.values,25))
    
    H_L_med=np.append(H_L_med,np.nanpercentile(Data.H.values/Data.L.values,50))
    H_L_top=np.append(H_L_top,np.nanpercentile(Data.H.values/Data.L.values,75))
    H_L_low=np.append(H_L_low,np.nanpercentile(Data.H.values/Data.L.values,25))
    
    H_l_med=np.append(H_l_med,np.nanpercentile(Data.H.values/Data.l.values,50))
    H_l_top=np.append(H_l_top,np.nanpercentile(Data.H.values/Data.l.values,75))
    H_l_low=np.append(H_l_low,np.nanpercentile(Data.H.values/Data.l.values,25))
    
    l_med=np.append(l_med,np.nanpercentile(Data.l.values,50))
    l_top=np.append(l_top,np.nanpercentile(Data.l.values,75))
    l_low=np.append(l_low,np.nanpercentile(Data.l.values,25))
    
reals=~np.isnan(H_all+L_all+l_all)
L_violation=np.sum((L_all[reals]<L_lim[0])+(L_all[reals]>L_lim[1]))/np.sum(reals)*100
H_L_violation=np.sum(H_all[reals]/L_all[reals]>H_L_max)/np.sum(reals)*100
H_l_violation=np.sum(H_all[reals]/l_all[reals]>H_l_max)/np.sum(reals)*100
l_pblh_violation=np.sum(l_all[reals]/pblh>l_pblh_max)/np.sum(reals)*100

print(f'L violated {L_violation}% of the times')
print(f'H/L violated {H_L_violation}% of the times')
print(f'H/l violated {H_l_violation}% of the times')
print(f'l/PBLH violated {l_pblh_violation}% of the times')

#%% Plots
#histogram
fig=plt.figure(figsize=(18,3))
ax=plt.subplot(1,4,1)
plt.hist(L_all[reals],25,color='k')
plt.plot([L_lim[0],L_lim[0]],[0,np.sum(reals)*0.75],'--r',linewidth=2)
plt.plot([L_lim[1],L_lim[1]],[0,np.sum(reals)*0.75],'--r',linewidth=2)
plt.xlabel(r'$L$ [m]')
plt.grid()

ax=plt.subplot(1,4,2)
plt.hist(H_all[reals]/L_all[reals],25,color='k')
plt.plot([H_L_max,H_L_max],[0,np.sum(reals)*0.75],'--r',linewidth=2)
plt.xlabel(r'$H/L$ [m]')
plt.grid()

ax=plt.subplot(1,4,3)
plt.hist(H_all[reals]/l_all[reals],25,color='k')
plt.plot([H_l_max,H_l_max],[0,np.sum(reals)*0.75],'--r',linewidth=2)
plt.xlabel(r'$H/l$ [m]')
plt.grid()

ax=plt.subplot(1,4,4)
plt.hist(l_all[reals]/pblh,25,color='k')
plt.plot([l_pblh_max,l_pblh_max],[0,np.sum(reals)*0.75],'--r',linewidth=2)
plt.xlabel(r'$l/PBLH$ [m]')
plt.grid()

utl.remove_labels(fig)
plt.tight_layout()

#polar plots
th=np.arange(0,2*np.pi,np.pi/100)
fig=plt.figure(figsize=(18,3))
ax=plt.subplot(1,4,1,polar=True)
ax.plot(wrap(wd/180*np.pi),wrap(L_med),color='k',marker='.',markersize=10)
ax.fill_between(wrap(wd/180*np.pi), wrap(L_low), wrap(L_top), color='k', alpha=0.5)
ax.plot(th,th*0+L_lim[0],'--r')
ax.plot(th,th*0+L_lim[1],'--r')
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi/2)
ax.set_rscale('symlog',linthresh=0.5)
ax.set_rlim([10,20000])
ax.set_rticks([],labels=[])
ax.set_thetagrids([0,45,90,135,180,225,270,315],labels=['N','NE','E','SE','S','SW','W','NW'])

ax=plt.subplot(1,4,2,polar=True)
ax.plot(wrap(wd/180*np.pi),wrap(H_L_med),color='k',marker='.',markersize=10)
ax.fill_between(wrap(wd/180*np.pi), wrap(H_L_low), wrap(H_L_top), color='k', alpha=0.5)
ax.plot(th,th*0+H_L_max,'--r')
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi/2)
ax.set_rlim([0,0.1])
ax.set_rticks([],labels=[])
ax.set_thetagrids([0,45,90,135,180,225,270,315],labels=['N','NE','E','SE','S','SW','W','NW'])

ax=plt.subplot(1,4,3,polar=True)
ax.plot(wrap(wd/180*np.pi),wrap(H_l_med),color='k',marker='.',markersize=10)
ax.fill_between(wrap(wd/180*np.pi), wrap(H_l_low), wrap(H_l_top), color='k', alpha=0.5)
ax.plot(th,th*0+H_l_max,'--r')
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi/2)
ax.set_rlim([0,2])
ax.set_rticks([],labels=[])
ax.set_thetagrids([0,45,90,135,180,225,270,315],labels=['N','NE','E','SE','S','SW','W','NW'])

ax=plt.subplot(1,4,4,polar=True)
ax.plot(wrap(wd/180*np.pi),wrap(l_med)/pblh,color='k',marker='.',markersize=10)
ax.fill_between(wrap(wd/180*np.pi), wrap(l_low)/pblh, wrap(l_top)/pblh, color='k', alpha=0.5)
ax.plot(th,th*0+l_pblh_max,'--r')
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi/2)
ax.set_rlim([0,0.3])
ax.set_rticks([],labels=[])
ax.set_thetagrids([0,45,90,135,180,225,270,315],labels=['N','NE','E','SE','S','SW','W','NW'])


