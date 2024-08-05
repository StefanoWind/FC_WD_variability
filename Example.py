# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 09:33:50 2024

@author: sletizia
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
import pandas as pd
import utm
from scipy.interpolate import griddata

warnings.filterwarnings('ignore')
plt.close('all')

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.size'] = 14

#%% Inputs
source=os.path.join(cd,'data/FC_topo_v2.nc')
source_nwtc='data/NWTC.xlsx'

wd=270#[deg] wind direction
y_sel=0#[m] selected profile

tapering=1000#[m] tapering length
response_cutoff=0.9 #value of smoothing filter response at cutoff
L_cutoff=1000#[m] smoothing cutoff lenghscale

#%% Initialization
Data=xr.open_dataset(source)

FC=pd.read_excel(source_nwtc).set_index('Site')

#locations
xy=utm.from_latlon(FC['Lat'].values,FC['Lon'].values)
xref=utm.from_latlon(FC['Lat']['M2'],FC['Lon']['M2'])[0]
yref=utm.from_latlon(FC['Lat']['M2'],FC['Lon']['M2'])[1]
FC['x']=xy[0]-xref
FC['y']=xy[1]-yref

x=Data.x.values-xref
y=Data.y.values-yref
X,Y=np.meshgrid(x,y)
Z=Data.z.values.T

sigma=(-np.log(response_cutoff)/2)**0.5/np.pi*L_cutoff#[m] smoothing sigma

#%% Main

#rotation
X_rot= X*utl.cosd(270-wd)+Y*utl.sind(270-wd)
Y_rot=-X*utl.sind(270-wd)+Y*utl.cosd(270-wd)

points=np.array([X_rot.ravel(),Y_rot.ravel()]).T
values=Z.ravel()
Z_rot=griddata(points,values,(X,Y))

#transformation
Z_trans_rot=Z_rot*0+np.nan
for iy in range(len(y)):
    z=Z_rot[iy,:]
    
    reals=~np.isnan(z)
    x_sel=x[reals]
    z_sel=z[reals]
    
    #detrending
    lf=np.polyfit(x_sel, z_sel,1)
    z_det=z_sel-lf[0]*x_sel-lf[1]
    
    #tapering matrix (if nans are present)
    T=np.eye(len(x_sel))
    for i in range(len(x_sel)):
        if x_sel[i]-x_sel[0]<tapering:
            T[i,i]=3*((x_sel[i]-x_sel[0])/tapering)**2-2*((x_sel[i]-x_sel[0])/tapering)**3
        if x_sel[-1]-x_sel[i]<tapering:
            T[i,i]=3*((x_sel[-1]-x_sel[i])/tapering)**2-2*((x_sel[-1]-x_sel[i])/tapering)**3

    #smoothing matrix (if nans are present)
    M=np.zeros((len(x_sel),len(x_sel)))
    for j in range(len(x_sel)):
        M[j,:]=np.exp(-(x_sel-x_sel[j])**2/(2*sigma**2))
        M[j,:]=M[j,:]/np.sum(M[j,:])

    z_smooth=np.matmul(M,z_det)#smoothing
    z_taper=np.matmul(T,z_smooth-np.nanmin(z_smooth))#tapering
    
    Z_trans_rot[iy,reals]=z_taper
    if np.abs(y[iy]-y_sel)<10:
        z_sel2=z+0.0
        z_det_sel=z*np.nan
        z_det_sel[reals]=z_det.copy()
        z_smooth_sel=z*np.nan
        z_smooth_sel[reals]=z_smooth.copy()
        z_taper_sel=z*np.nan
        z_taper_sel[reals]=z_taper.copy()
        
    print(iy/len(y))

X_antirot= X*utl.cosd(270-wd)-Y*utl.sind(270-wd)
Y_antirot= X*utl.sind(270-wd)+Y*utl.cosd(270-wd)

points=np.array([X_antirot.ravel(),Y_antirot.ravel()]).T
values=Z_trans_rot.ravel()
Z_trans=griddata(points,values,(X,Y))
 
 
#%% Main
plt.close('all')
fig=plt.figure(figsize=(10,10))
ax=fig.add_subplot(111)
plt.pcolor(Data.x-xref,Data.y-yref,Data.z.T,cmap='gray',vmin=1800,vmax=1880)
plt.xlim([-100,1900])
plt.ylim([-1000,1000])
utl.axis_equal()

fig=plt.figure(figsize=(18,6))
ax=fig.add_subplot(111)
ax.fill_between(x,x*0,z_sel2,color='k',alpha=0.25)
ax.fill_between(x,x*0,z_det_sel,color='k',alpha=1)

fig=plt.figure(figsize=(18,6))
ax=fig.add_subplot(111)
ax.fill_between(x,x*0,z_det_sel,color='k',alpha=0.25)
ax.fill_between(x,x*0,z_smooth_sel,color='k',alpha=1)

fig=plt.figure(figsize=(18,6))
ax=fig.add_subplot(111)
ax.fill_between(x,x*0,z_smooth_sel,color='k',alpha=0.25)
ax.fill_between(x,x*0,z_taper_sel,color='k',alpha=1)

fig=plt.figure(figsize=(10,10))
ax=fig.add_subplot(111)
plt.pcolor(Data.x-xref,Data.y-yref,Z_trans,cmap='gray',vmin=50,vmax=110)
plt.xlim([-100,1900])
plt.ylim([-1000,1000])
utl.axis_equal()