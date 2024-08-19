# -*- coding: utf-8 -*-
"""
Compare speed up with field data
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
import glob

warnings.filterwarnings('ignore')
plt.close('all')

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.size'] = 10

#%% Inputs
source=os.path.join(cd,'data/20deg/DS*nc')
source_nwtc='data/NWTC.xlsx'
sources_ws={'low':'data/wspdDifference_wdirBins_M2 5m_M5 3m_20230101_20240406.csv',
            'high':'data/wspdDifference_wdirBins_M2 50m_M5 61m_20230401_20240508.csv'}
site_ref='M2'
site_sel='M5'

z_ref={'low':5.0,
       'high': 50.0}
z_sel={'low':3.0,
       'high': 61.0}

z0=0.03#[m] roughness length

#graphics 
xlim=[-250,1500]#[m]
ylim=[-875,875]#[m]
colors=['r','b','g']
labels={'low':    'M2 5m and M5 3m',
        'high':'M2 50m and M5 61m'}


#%% Initalization
FC=pd.read_excel(source_nwtc).set_index('Site')

#locations
xy=utm.from_latlon(FC['Lat'].values,FC['Lon'].values)
xref=utm.from_latlon(FC['Lat']['M2'],FC['Lon']['M2'])[0]
yref=utm.from_latlon(FC['Lat']['M2'],FC['Lon']['M2'])[1]
FC['x']=xy[0]-xref
FC['y']=xy[1]-yref

#graphics
files=glob.glob(source)
fig=plt.figure(figsize=(18,10))
gs = gridspec.GridSpec(2, int(np.ceil(len(files)/2))+1, width_ratios=[1]*int(np.ceil(len(files)/2))+[0.1])

#zeroing
wd_all=[]
DS_ref_all=[]
DS_sel_all=[]

#%% Main

for f in files:
    Data=xr.open_dataset(f)
    
    wd=np.float64(f[-6:-3])

    DS_ref=Data.DS.interp({'x':FC['x'][site_ref],'y':FC['y'][site_ref]})/100
    DS_sel=Data.DS.interp({'x':FC['x'][site_sel],'y':FC['y'][site_sel]})/100
    # print('Coordinate swapped')
    DS_ref_all=np.append(DS_ref_all,DS_ref)
    DS_sel_all=np.append(DS_sel_all,DS_sel)
    wd_all=np.append(wd_all,wd)



#%% Plots
matplotlib.rcParams['font.size'] = 22
plt.figure(figsize=(18,8))
ctr=0
for s in sources_ws.keys():
    DWS=((1+DS_sel_all)/(1+DS_ref_all)*np.log(z_sel[s]/z0)/np.log(z_ref[s]/z0)-1)*100
   
    WS=pd.read_csv(sources_ws[s])
    print(s+': corr = ' +str(utl.nancorrcoef(WS['Wind speed difference (m/s)'],DWS)[0,1]))
    plt.plot(WS['Wind direction (deg)'],WS['Wind speed difference (m/s)'],'.-',color=colors[ctr],markersize=10,label='Exp. ('+labels[s]+')')
    plt.plot(wd_all,DWS,'--',color=colors[ctr],markersize=10,label='JH75 ('+labels[s]+')')
    ctr+=1
DWS_noshear=((1+DS_sel_all)/(1+DS_ref_all)-1)*100
plt.plot(wd_all,DWS_noshear,'--k',markersize=10,label='JH75 (no shear)')
plt.legend(draggable=True)
plt.xlabel('Wind direction (deg)')
plt.ylabel('Relative wind speed difference (percent)')
plt.xticks(wd_all)
plt.grid()
plt.tight_layout()


# matplotlib.rcParams['font.size'] = 22
# plt.figure(figsize=(18,6))
# ctr=0

# WS=pd.read_csv(sources_ws[1])
# rho=utl.nancorrcoef(WS['Wind speed difference (m/s)'],DWS)[0,1]
# plt.plot(WS['Wind direction (deg)'],WS['Wind speed difference (m/s)']-WS['Wind speed difference (m/s)'].mean(),'.-',color=colors[ctr],markersize=10,label=labels[sources_ws[1]]+r', $\rho='+str(np.round(rho,2))+'$')
# plt.plot(wd_all,DWS,'.-k',markersize=10,label='Linear momentum theory')
# plt.legend()
# plt.xlabel(r'$\theta_w$ [$^\circ$]')
# plt.ylabel(r'$\frac{\overline{u}_{M2}-\overline{u}_{M5}}{\overline{u}_{M5}}$ [%]')
# plt.xticks(wd_all)
# plt.grid()
# plt.tight_layout()
