# -*- coding: utf-8 -*-
"""
Plot NWTC layout
"""
import os
cd=os.path.dirname(__file__)
import sys
sys.path.append('C:/Users/SLETIZIA/OneDrive - NREL/Desktop/PostDoc/utils')
import utils as utl
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import matplotlib
import rasterio
import utm
import xarray as xr
from scipy.interpolate import griddata
from matplotlib.patches import Polygon

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.size'] = 14
plt.close('all')

#%% Inputs
source_topo='data/USGS_13_n40w106_20230602.tif'
source_nwtc='data/NWTC.xlsx'

margin=0.125#[deg] margin around locations in degrees
skip=1#skip topography points
dx_topo=10#[m] W-E resolution
dy_topo=10#[m] S-N resolution
RIX_slope=17#[deg]

#graphics
xlim1=[-10000,10000]#[m]
ylim1=[-10000,10000]#[m]

xlim2=[-250,1500]#[m]
ylim2=[-875,875]#[m]

xlim3=[-5000,5000]#[m]
ylim3=[-5000,5000]#[m]

skip1=int(500/dx_topo)
skip2=int(100/dx_topo)

#%% Initialization
file = rasterio.open(source_topo)
NWTC=pd.read_excel(source_nwtc).set_index('Site')

#%% Main

#locations
xy=utm.from_latlon(NWTC['Lat'].values,NWTC['Lon'].values)
xref=xy[0][0]
yref=xy[1][0]
NWTC['x']=xy[0]-xref
NWTC['y']=xy[1]-yref

#read topography
print('Loading topography')
Z_topo0 = file.read()[0][::skip,::skip]
cols, rows = np.meshgrid(np.arange(0,file.width,skip), np.arange(0,file.height,skip))
LON0, LAT0 = rasterio.transform.xy(file.transform, rows, cols)
sel=(np.array(LON0)>np.min(NWTC['Lon'])-margin)*(np.array(LON0)<np.max(NWTC['Lon'])+margin)*(np.array(LAT0)>np.min(NWTC['Lat'])-margin)*(np.array(LAT0)<np.max(NWTC['Lat'])+margin)

#interpolate topography on regular grid
print('Interpolating topography')
XY=utm.from_latlon(np.array(LAT0),np.array(LON0))
X_topo0=XY[0]
Y_topo0=XY[1]

points=np.zeros((np.sum(sel),2))
points[:,0]=X_topo0[sel]
points[:,1]=Y_topo0[sel]
values=Z_topo0[sel]
x_topo=np.arange(np.min(X_topo0[sel]),np.max(X_topo0[sel]),dx_topo)
y_topo=np.arange(np.min(Y_topo0[sel]),np.max(Y_topo0[sel]),dy_topo)
X_topo,Y_topo=np.meshgrid(x_topo,y_topo)
Z_topo=griddata(points,values,(X_topo,Y_topo))

#altitude
NWTC['z']=griddata(points,values,(NWTC['x']+xref,NWTC['y']+yref))
zref=NWTC['z']['M2']

#slope
Topo=xr.DataArray(data=Z_topo.T,coords={'x':x_topo,'y':y_topo})
dZdx= Topo.differentiate(coord='x').values.T
dZdy= Topo.differentiate(coord='y').values.T
slope=utl.arctand((dZdx**2+dZdy**2)**0.5)
direction=utl.cart2pol(dZdx, dZdy)[1]

#solar array location
x_solar=0
y_solar=0
N=0
vertices=[]
for s in NWTC.index:
    if 'Solar field' in s:
        x_solar+=NWTC['x'][s]
        y_solar+=NWTC['y'][s]
        N+=1
        vertices.append((NWTC['x'][s],NWTC['y'][s]))

x_solar=x_solar/N
y_solar=y_solar/N

#turbine locations
x_turbine=[]
y_turbine=[]
for s in NWTC.index:
    if 'turbine' in s:
        x_turbine.append(NWTC['x'][s])
        y_turbine.append(NWTC['y'][s])

#distances
print('M2-M5: '+str(np.round(((NWTC['x']['M2']-NWTC['x']['M5'])**2+(NWTC['y']['M2']-NWTC['y']['M5'])**2)**0.5)))
print('M2-Solar field: '+str(np.round(((NWTC['x']['M2']-x_solar)**2+(NWTC['y']['M2']-y_solar)**2)**0.5)))
print('M5-Solar field: '+str(np.round(((NWTC['x']['M5']-x_solar)**2+(NWTC['y']['M5']-y_solar)**2)**0.5)))

sel_x=(x_topo-xref>xlim2[0])*(x_topo-xref<xlim2[1])
sel_y=(y_topo-yref>ylim2[0])*(y_topo-yref<ylim2[1])
print('Median slope ='+str(np.nanmedian(slope[sel_y,:][:,sel_x])))

RIX=np.sum(slope[sel_y,:][:,sel_x]>RIX_slope)/np.sum(~np.isnan(slope[sel_y,:][:,sel_x]))*100
print('RIX = '+str(RIX))

#%% Output
sel_x=(x_topo-xref>xlim3[0])*(x_topo-xref<xlim3[1])
sel_y=(y_topo-yref>ylim3[0])*(y_topo-yref<ylim3[1])
Output=xr.Dataset()
Output['z']=xr.DataArray(data=Z_topo[sel_y,:][:,sel_x].T,coords={'x':x_topo[sel_x],'y':y_topo[sel_y]})
Output.to_netcdf(os.path.join(cd,'data/FC_topo_v2.nc'))

#%% Plots
print('Plotting')
plt.close('all')
plt.figure(figsize=(14,10))
ax=plt.subplot(2,2,1)
plt.contourf(x_topo-xref,y_topo-yref,Z_topo,np.arange(1600,2501,20),cmap='summer',vmin=1600,vmax=2500,extend='both')
polygon = Polygon([(xlim2[0],ylim2[0]),(xlim2[0],ylim2[1]),(xlim2[1],ylim2[1]),(xlim2[1],ylim2[0])], closed=True, edgecolor='k',linewidth=2,facecolor=(0,0,0,0))
ax.add_patch(polygon)
plt.xlim(xlim1)
plt.ylim(ylim1)
plt.xticks(np.arange(xlim1[0],xlim1[1]+1,5000))
plt.yticks(np.arange(ylim1[0],ylim1[1]+1,5000))
plt.xlabel('W-E [m]')
plt.ylabel('S-N [m]')
plt.grid()
plt.colorbar(label='Altitude above sea level [m]',ticks=np.arange(1600,2501,100))
utl.axis_equal()

ax=plt.subplot(2,2,2)
plt.contourf(x_topo-xref,y_topo-yref,slope,np.arange(0,31,2.5),cmap='summer',vmin=0,vmax=30,extend='both')
polygon = Polygon([(xlim2[0],ylim2[0]),(xlim2[0],ylim2[1]),(xlim2[1],ylim2[1]),(xlim2[1],ylim2[0])], closed=True, edgecolor='k',linewidth=2,facecolor=(0,0,0,0))
ax.add_patch(polygon)
plt.xlim(xlim1)
plt.ylim(ylim1)
plt.xticks(np.arange(xlim1[0],xlim1[1]+1,5000))
plt.yticks(np.arange(ylim1[0],ylim1[1]+1,5000))
plt.xlabel('W-E [m]')
plt.ylabel('S-N [m]')
plt.grid()
plt.colorbar(label='Slope magnitude [$^\circ$]',ticks=np.arange(0,31,3))
plt.quiver(X_topo[::skip1,::skip1]-xref,Y_topo[::skip1,::skip1]-yref,utl.cosd(direction[::skip1,::skip1]),utl.sind(direction[::skip1,::skip1]),scale=50,color=(0,0,0,0.5))
utl.axis_equal()

ax=plt.subplot(2,2,3)
plt.contourf(x_topo-xref,y_topo-yref,Z_topo,np.arange(1810,1866),cmap='summer',vmin=1820,vmax=1865,extend='both')
plt.plot(NWTC['x']['M2'],NWTC['y']['M2'],'^k',markersize=15)
plt.plot(NWTC['x']['M5'],NWTC['y']['M5'],'^k',markersize=15)
plt.plot(x_turbine,y_turbine,'ow',markersize=15,linewidth=5)
plt.plot(x_turbine,y_turbine,'3k',markersize=15)
plt.plot(x_turbine,y_turbine,'3k',markersize=15,linewidth=2)
polygon = Polygon(vertices, closed=True, edgecolor='k',facecolor=(0,0,0,0.5))
ax.add_patch(polygon)
plt.plot(x_solar,y_solar,'xk',markersize=5)
plt.plot([NWTC['x']['M2'],NWTC['x']['M5']],[NWTC['y']['M2'],NWTC['y']['M5']],'--k')
plt.plot([NWTC['x']['M2'],x_solar],[NWTC['y']['M2'],y_solar],'--k')
plt.plot([NWTC['x']['M5'],x_solar],[NWTC['y']['M5'],y_solar],'--k')          
plt.xlim(xlim2)
plt.ylim(ylim2)
plt.xticks(np.arange(xlim2[0],xlim2[1]+1,250))
plt.yticks(np.arange(ylim2[0],ylim2[1]+1,250))
plt.xlabel('W-E [m]')
plt.ylabel('S-N [m]')
plt.grid()
plt.colorbar(label='Altitude above sea level [m]',ticks=np.arange(1810,1866,5))
utl.axis_equal()

ax=plt.subplot(2,2,4)
plt.contourf(x_topo-xref,y_topo-yref,slope,np.arange(0,31,2.5),cmap='summer',vmin=0,vmax=30,extend='both')
plt.plot(NWTC['x']['M2'],NWTC['y']['M2'],'^k',markersize=15)
plt.plot(NWTC['x']['M5'],NWTC['y']['M5'],'^k',markersize=15)
plt.plot(x_turbine,y_turbine,'ow',markersize=15,linewidth=5)
plt.plot(x_turbine,y_turbine,'3k',markersize=15)
polygon = Polygon(vertices, closed=True, edgecolor='k',facecolor=(0,0,0,0.5))
ax.add_patch(polygon)
plt.xlim(xlim2)
plt.ylim(ylim2)
plt.xticks(np.arange(xlim2[0],xlim2[1]+1,250))
plt.yticks(np.arange(ylim2[0],ylim2[1]+1,250))
plt.xlabel('W-E [m]')
plt.ylabel('S-N [m]')
plt.grid()
plt.colorbar(label='Slope magnitude [$^\circ$]',ticks=np.arange(0,31,3))
plt.quiver(X_topo[::skip2,::skip2]-xref,Y_topo[::skip2,::skip2]-yref,utl.cosd(direction[::skip2,::skip2]),utl.sind(direction[::skip2,::skip2]),scale=25,color=(0,0,0,0.5))
utl.axis_equal()