# -*- coding: utf-8 -*-
"""
Plots maps do speedup for different wind directions
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
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon
import glob

warnings.filterwarnings('ignore')
plt.close('all')

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.size'] = 14

#%% Inputs
source=os.path.join(cd,'data/DS*nc')
source_nwtc='data/NWTC.xlsx'
source_ws='data/wspdDifference_wdirBins_M5 61m_M2 50m_20230401_20240508.csv'

#graphics 
xlim=[-250,1500]#[m]
ylim=[-875,875]#[m]

#%% Initalization
FC=pd.read_excel(source_nwtc).set_index('Site')

WS=pd.read_csv(source_ws)

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
x_solar=0
y_solar=0
N=0
vertices=[]
for s in FC.index:
    if 'Solar field' in s:
        x_solar+=FC['x'][s]
        y_solar+=FC['y'][s]
        N+=1
        vertices.append((FC['x'][s],FC['y'][s]))

x_solar=x_solar/N
y_solar=y_solar/N

files=glob.glob(source)
fig=plt.figure(figsize=(18,10))
gs = gridspec.GridSpec(2, int(np.ceil(len(files)/2))+1, width_ratios=[1]*int(np.ceil(len(files)/2))+[0.1])

wd_all=[]
DS_selected_all=[]

#%% Main
ctr=0
for f in files:
    Data=xr.open_dataset(f)
    wd=np.float64(f[-6:-3])
        
    sel_x=(Data.x>xlim[0])*(Data.x<xlim[1])
    sel_y=(Data.y>ylim[0])*(Data.y<ylim[1])
    
    X,Y=np.meshgrid(Data.x[sel_x],Data.y[sel_y])
    
    points=np.array([X.ravel(),Y.ravel()]).T
    values=Data.DS.values[sel_y,:][:,sel_x].ravel()
    DS_selected=griddata(points,values,([FC['x']['M2'],FC['x']['M5'],x_solar],[FC['y']['M2'],FC['y']['M5'],y_solar]))
    
    ax = plt.subplot(gs[int(ctr/int(np.ceil(len(files)/2))), ctr-int(ctr/int(np.ceil(len(files)/2)))*int(np.ceil(len(files)/2))])
   
    cf=plt.contourf(Data.x,Data.y,Data.DS,np.arange(-30,31,2.5),vmin=-30,vmax=30,cmap='seismic',extend='both')
    
    plt.plot(FC['x']['M2'],FC['y']['M2'],'^g',markersize=15)
    plt.text(FC['x']['M2']-150,FC['y']['M2']-250,str(np.round(DS_selected[0],1))+' %',bbox=dict(facecolor='white', alpha=0.5))
    
    plt.plot(FC['x']['M5'],FC['y']['M5'],'^g',markersize=15)
    plt.text(FC['x']['M5']-150,FC['y']['M5']-250,str(np.round(DS_selected[1],1))+' %',bbox=dict(facecolor='white', alpha=0.5))
    
   
    polygon = Polygon(vertices, closed=True, edgecolor='g',facecolor=(0,1,0,0.5))
    ax.add_patch(polygon)     
    plt.text(x_solar-150,y_solar-250,str(np.round(DS_selected[2],1))+' %',bbox=dict(facecolor='white', alpha=0.5))
    
    plt.plot(x_turbine,y_turbine,'ow',markersize=15,linewidth=5)
    plt.plot(x_turbine,y_turbine,'3k',markersize=15)
    plt.plot(x_turbine,y_turbine,'3k',markersize=15,linewidth=2)
    
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xticks(np.arange(xlim[0],xlim[1]+1,250),rotation=30)
    plt.yticks(np.arange(ylim[0],ylim[1]+1,250))
    plt.xlabel('W-E [m]')
    plt.ylabel('S-N [m]')
    plt.grid()
    ax.arrow(1250,625,utl.cosd(270-wd)*150,utl.sind(270-wd)*150,head_width=300, head_length=100, fc='g', ec='k',width=200)
    utl.axis_equal()
    ctr+=1
    
    DS_selected_all=utl.vstack(DS_selected_all,DS_selected)
    wd_all=np.append(wd_all,wd)

utl.remove_labels(fig)

cbar_ax = plt.subplot(gs[:, -1])
cbar = fig.colorbar(cf, cax=cbar_ax,label=r'$\Delta S$ [%]')

#%% Plots
matplotlib.rcParams['font.size'] = 22
plt.figure(figsize=(18,6))
plt.plot(WS['Wind direction (deg)'],WS['Wind speed difference (m/s)'],'.-b',markersize=10,label='Data')
plt.plot(wd_all,(DS_selected_all[:,0]/100-DS_selected_all[:,1]/100)/(DS_selected_all[:,1]/100+1)*100,'.-k',markersize=10,label='JH75 model')
plt.legend()
plt.xlabel(r'$\theta_w$ [$^\circ$]')
plt.ylabel(r'$\frac{\overline{u}_{M2}(z=h+l)-\overline{u}_{M5}(z=h+l)}{\overline{u}_{M5}(z=h+l)}$ [%]')
plt.xticks(wd_all)
plt.grid()
plt.tight_layout()

