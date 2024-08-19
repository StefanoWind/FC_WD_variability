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
source=os.path.join(cd,'data/45deg/DS*nc')
source_nwtc='data/NWTC.xlsx'

#graphics 
xlim=[-250,1500]#[m]
ylim=[-875,875]#[m]
plot_rows=2
plot_cols=4

#%% Initalization
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
gs = gridspec.GridSpec(plot_rows, plot_cols+1, width_ratios=[1]*plot_cols+[0.1])

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
    values=Data.DS.values[sel_x,:][:,sel_y].T.ravel()
    DS_sel2=griddata(points,values,([FC['x']['M2'],FC['x']['M5'],x_solar],[FC['y']['M2'],FC['y']['M5'],y_solar]))
    
    DS_sel=[np.float64(Data.DS.interp({'x':FC['x']['M2'],'y':FC['y']['M2']}).values),
            np.float(Data.DS.interp({'x':FC['x']['M5'],'y':FC['y']['M5']}).values),
            np.float(Data.DS.interp({'x':x_solar,'y':y_solar}).values)]
    
    print(str(wd)+': '+str(np.round(DS_sel,2)))
    print(str(wd)+': '+str(np.round(DS_sel2,2)))
    
    ax = plt.subplot(gs[int(ctr/plot_cols), ctr-int(ctr/plot_cols)*plot_cols])
   
    cf=plt.contourf(Data.x,Data.y,Data.DS.T,np.arange(-30,31,2.5),vmin=-30,vmax=30,cmap='seismic',extend='both')
    
    plt.plot(FC['x']['M2'],FC['y']['M2'],'^g',markersize=15)
    plt.text(FC['x']['M2']-150,FC['y']['M2']-250,str(np.round(DS_sel[0],1))+' %',bbox=dict(facecolor='white', alpha=0.5))
    
    plt.plot(FC['x']['M5'],FC['y']['M5'],'^g',markersize=15)
    plt.text(FC['x']['M5']-150,FC['y']['M5']-250,str(np.round(DS_sel[1],1))+' %',bbox=dict(facecolor='white', alpha=0.5))
   
    polygon = Polygon(vertices, closed=True, edgecolor='g',facecolor=(0,1,0,0.5))
    ax.add_patch(polygon)     
    plt.text(x_solar-150,y_solar-250,str(np.round(DS_sel[2],1))+' %',bbox=dict(facecolor='white', alpha=0.5))
    
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

utl.remove_labels(fig)

cbar_ax = plt.subplot(gs[:, -1])
cbar = fig.colorbar(cf, cax=cbar_ax,label=r'$\Delta S$ [%]')



