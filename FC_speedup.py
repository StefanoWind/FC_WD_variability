# -*- coding: utf-8 -*-
'''
Hill speed up from Jackosn and Hunt 1975 based on real FC topography
'''
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
import terrain_utils as trn
import utm
from scipy.interpolate import griddata
from matplotlib.patches import Polygon

warnings.filterwarnings('ignore')
plt.close('all')

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.size'] = 14

#%% Inputs
source=os.path.join(cd,'data/FC_topo_v2.nc')
source_nwtc='data/NWTC.xlsx'

z0=0.03#[m] surface roughness (farmland with very few buildings/trees in WAsP)
wds=np.arange(10,360,20)#[deg] set of wind directions

xlim=[-250,1500]#[m]
ylim=[-875,875]#[m]

tapering=1000#[m] tapering length
response_cutoff=0.9 #value of smoothing filter response at cutoff
L_cutoff=100#[m] smoothing cutoff lenghscale

#graphics
dist_plot=200#[m]

#%% Initalization
Data=xr.open_dataset(source)
FC=pd.read_excel(source_nwtc).set_index('Site')

#locations
xy=utm.from_latlon(FC['Lat'].values,FC['Lon'].values)
xref=utm.from_latlon(FC['Lat']['M2'],FC['Lon']['M2'])[0]
yref=utm.from_latlon(FC['Lat']['M2'],FC['Lon']['M2'])[1]
FC['x']=xy[0]-xref
FC['y']=xy[1]-yref

#turbine locations
x_turbine=[]
y_turbine=[]
for s in FC.index:
    if 'turbine' in s:
        x_turbine.append(FC['x'][s])
        y_turbine.append(FC['y'][s])

#solar array location
vertices=[]
for s in FC.index:
    if 'Solar field' in s:
        vertices.append((FC['x'][s],FC['y'][s]))

#grid
x=Data.x.values-xref
y=Data.y.values-yref
Z=Data.z.values.T

dx=x[1]-x[0]
dy=y[1]-y[0]
X,Y=np.meshgrid(x,y)

sigma=(-np.log(response_cutoff)/2)**0.5/np.pi*L_cutoff#[m] smoothing sigma

#%% Main
for wd in wds:
    #rotation
    X_rot= X*utl.cosd(270-wd)+Y*utl.sind(270-wd)
    Y_rot=-X*utl.sind(270-wd)+Y*utl.cosd(270-wd)
    
    points=np.array([X_rot.ravel(),Y_rot.ravel()]).T
    values=Z.ravel()
    Z_rot=griddata(points,values,(X,Y))
    
    #speedup
    DS_rot=Z*0
    L_rot=Z*0
    H_rot=Z*0
    l_rot=Z*0
    for j in range(len(y)):
        reals=~np.isnan(Z_rot[j,:])
        if np.sum(reals)>1:
            DS_rot[j,:],f_all,H_rot[j,:],L_rot[j,:],l_rot[j,:]=trn.speedup(x,Z_rot[j,:],tapering,sigma,z0)
        else:
            DS_rot[j,:]=np.nan
        print(j/len(y))
        
    #sample profiles
    jref=np.argmin(np.abs(y))
    jsel=[jref-int(dist_plot/dy),jref,jref+int(dist_plot/dy)]
    
    #anti-rotation
    X_antirot= X*utl.cosd(270-wd)-Y*utl.sind(270-wd)
    Y_antirot= X*utl.sind(270-wd)+Y*utl.cosd(270-wd)
    
    points=np.array([X_antirot.ravel(),Y_antirot.ravel()]).T
    values=DS_rot.ravel()
    DS=griddata(points,values,(X,Y))
    
    values=H_rot.ravel()
    H_all=griddata(points,values,(X,Y))
   
    values=L_rot.ravel()
    L_all=griddata(points,values,(X,Y))
    
    values=l_rot.ravel()
    l_all=griddata(points,values,(X,Y))
    
    #Output
    Output=xr.Dataset()
    Output['DS']=xr.DataArray(data=DS.T,  coords={'x':x,'y':y})
    Output['L']=xr.DataArray(data=L_all.T,coords={'x':x,'y':y})
    Output['H']=xr.DataArray(data=H_all.T,coords={'x':x,'y':y})
    Output['l']=xr.DataArray(data=l_all.T,coords={'x':x,'y':y})
    Output.to_netcdf(os.path.join(cd,f'data/DS_{wd:03d}.nc'))

    #Plots     
    #check on rotation, boundary effects
    plt.figure(figsize=(18,6))
    ax=plt.subplot(1,2,1)
    plt.contourf(X,Y,Z_rot,np.arange(1500,2501),cmap='summer',vmin=1500,vmax=2500,extend='both')
    plt.xlabel('$x$ [m]')
    plt.ylabel('$y$ [m]')
    plt.colorbar(label=r'$z$ [m]')
    
    utl.axis_equal()
    ax=plt.subplot(1,2,2)
    plt.contourf(X,Y,DS_rot,np.arange(-50,51,5),vmin=-50,vmax=50,cmap='seismic',extend='both')
    plt.xlabel('$x$ [m]')
    plt.ylabel('$y$ [m]')
    plt.colorbar(label=r'$\Delta S$ [%]')
    utl.axis_equal()
    plt.tight_layout()
    
    #speedup ratio
    plt.figure()
    plt.contourf(X,Y,DS,np.arange(-20,21,2.5),vmin=-20,vmax=20,cmap='seismic',extend='both')
    plt.plot(FC['x']['M2'],FC['y']['M2'],'^g',markersize=15)
    plt.plot(FC['x']['M5'],FC['y']['M5'],'^g',markersize=15)
    plt.plot(x_turbine,y_turbine,'ow',markersize=15,linewidth=5)
    plt.plot(x_turbine,y_turbine,'3k',markersize=15)
    plt.plot(x_turbine,y_turbine,'3k',markersize=15,linewidth=2)
    for j in jsel:
        plt.plot(X_antirot[j,:],Y_antirot[j,:],'--k',linewidth=2)
    polygon = Polygon(vertices, closed=True, edgecolor='g',facecolor=(0,1,0,0.5))
    plt.gca().add_patch(polygon)     
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xticks(np.arange(xlim[0],xlim[1]+1,250))
    plt.yticks(np.arange(ylim[0],ylim[1]+1,250))
    plt.xlabel('W-E [m]')
    plt.ylabel('S-N [m]')
    plt.grid()
    plt.colorbar(label=r'$\Delta S$ [%]')
    utl.axis_equal()
    
    #selected profiles
    plt.figure(figsize=(18,8))
    ctr=1
    for j in jsel:
        ax=plt.subplot(len(jsel),1, ctr)
        DS2,f,H,L,l=trn.speedup(x,Z_rot[j,:],tapering,sigma,z0)
        ax.fill_between(x,x*0,Z_rot[j,:]-np.min(Z_rot[j,:]),color='k',alpha=0.25)
        ax.fill_between(x,x*0,f*H,color='b',alpha=0.25)
        plt.plot(x,x*0+500,'-r',linewidth=1)
        ax.fill_between(x,x*0+500-20,x*0+500+20,color='r',alpha=0.25)
        plt.plot(x,DS2*2+500,'r')
        plt.xlabel(r'$x$ [m]')
        plt.ylabel(r'$z$ [m]')
        plt.xlim([x[0],x[-1]])
        plt.ylim([0,800])
        plt.text(x[0]+100,700,'ABCDEFG'[ctr-1])
        plt.grid()
        ctr+=1
    utl.save_all_fig(f'FC.speedup.{wd:03d}.', cd)
    plt.close('all')