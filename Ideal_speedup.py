# -*- coding: utf-8 -*-
'''
Hill speed up from Jackosn and Hunt 1975 based on real FC topography
'''
#2024-07-22: created, finalized
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
import terrain_utils as trn
from scipy.interpolate import griddata
warnings.filterwarnings('ignore')
plt.close('all')

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.size'] = 16

#%% Inputs
source=os.path.join(cd,'data/Ideal_hill.nc')

z0=0.03#[m] surface roughness (farmland with very few buildings/trees in WAsP)
wd=270-45#[deg]

tapering=1000#[m] tapering length
response_cutoff=0.9 #value of smoothing filter response at cutoff
L_cutoff=100#[m] smoothing cutoff lenghscale

#%% Initalization
Data=xr.open_dataset(source)
x=Data.x.values
y=Data.y.values
Z=Data.z.values.T

dx=x[1]-x[0]
dy=y[1]-y[0]
X,Y=np.meshgrid(x,y)

l0=1/8*z0*(np.float64(Data.L)/z0)**0.9
DS_max=np.float64(Data.H)/np.float64(Data.L)*np.log(np.float64(Data.L)/z0)**2/np.log(l0/z0)**2*100
DS_theory=np.float64(Data.H)/np.float64(Data.L)*np.log(np.float64(Data.L)/z0)**2/np.log(l0/z0)**2*(1-(x/np.float64(Data.L))**2)/(1+(x/np.float64(Data.L))**2)**2*100

sigma=(-np.log(response_cutoff)/2)**0.5/np.pi*L_cutoff#[m]

#%% Main

#rotation
X_rot= X*utl.cosd(270-wd)+Y*utl.sind(270-wd)
Y_rot=-X*utl.sind(270-wd)+Y*utl.cosd(270-wd)

points=np.array([X_rot.ravel(),Y_rot.ravel()]).T
values=Z.ravel()
Z_rot=griddata(points,values,(X,Y))

DS_rot=Z*0
f_rot=Z*0
H=np.zeros(len(y))
L=np.zeros(len(y))
l=np.zeros(len(y))
for j in range(len(y)):
    reals=~np.isnan(Z_rot[j,:])
    if np.sum(reals)>1:
        DS_rot[j,:],f_rot[j,:],H[j],L[j],l[j]=trn.speedup(x,Z_rot[j,:],tapering,sigma,z0)
    else:
        DS_rot[j,:]=np.nan
    print(j/len(y))
    
#anti-rotation
X_antirot= X*utl.cosd(270-wd)-Y*utl.sind(270-wd)
Y_antirot= X*utl.sind(270-wd)+Y*utl.cosd(270-wd)

points=np.array([X_antirot.ravel(),Y_antirot.ravel()]).T
values=DS_rot.ravel()
DS=griddata(points,values,(X,Y))

values=f_rot.ravel()
f=griddata(points,values,(X,Y))
    
#%% Plots
plt.close('all')
plt.figure(figsize=(18,4))
plt.subplot(1,2,1)
plt.contourf(X,Y,DS,np.arange(-DS_max,DS_max+0.1,0.1),vmin=-DS_max,vmax=DS_max,cmap='seismic',extend='both')
plt.xlabel('W-E [m]')
plt.ylabel('S-N [m]')
plt.grid()
plt.colorbar(label=r'$\Delta S$ [%]')
plt.title(r'$\Delta S_{max}='+str(np.round(DS_max,2))+' $ (theory), $\Delta S_{max}='+str(np.round(np.nanmax(DS),2))+'$')
utl.axis_equal()

plt.subplot(1,2,2)
plt.contourf(X,Y,f,np.arange(0,1.1,0.1),vmin=0,vmax=1,cmap='hot',extend='both')
plt.xlabel('W-E [m]')
plt.ylabel('S-N [m]')
plt.grid()
plt.colorbar(label=r'$f$')
utl.axis_equal()
plt.tight_layout()

plt.figure(figsize=(14,10))
plt.subplot(2,1,1)
plt.plot(x,DS_theory,'k',label='Hunt et al. 1988')
plt.plot(x,DS[y==0,:].squeeze(),'r',label='Numerical ['+str(np.round((np.nanmax(DS)/DS_max-1)*100,2))+'% max error]')
plt.legend()
plt.grid()
plt.xlabel(r'$x$ [m]')
plt.ylabel(r'$\Delta S$ [%]')
plt.title(r'$L='+str(np.float64(Data.L))+'$ m, $H='+str(np.float64(Data.H))+'$ m, $\sigma='+str(np.round(sigma,2))+'$ m')

plt.subplot(2,1,2)
plt.plot(x,Z[y==0,:].squeeze(),'k',label='Original')
plt.plot(x,f[y==0,:].squeeze()*H[y==0],'r',label='Rotated, tapered, smoothed')
plt.legend()
plt.grid()
plt.xlabel(r'$x$ [m]')
plt.ylabel(r'$z$ [m]')
plt.title(r'$L='+str(np.float64(Data.L))+'$ m, $H='+str(np.float64(Data.H))+'$ m, $\sigma='+str(np.round(sigma,2))+r'$ m, $\theta_w='+str(wd)+'^\circ$')
plt.tight_layout()
