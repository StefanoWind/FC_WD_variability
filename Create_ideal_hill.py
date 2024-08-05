# -*- coding: utf-8 -*-
"""
Create test topography
"""

import os
cd=os.getcwd()
import sys
sys.path.append('C:/Users/SLETIZIA/OneDrive - NREL/Desktop/PostDoc/utils')
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
import warnings
import matplotlib


warnings.filterwarnings('ignore')
plt.close('all')

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.size'] = 16


#%% Inputs
L=100#[m] length scales
H=5#[m] height
z0=0.03#[m] surface roughness (farmland with very few buildings/trees in WAsP)

#domain
xmin=-2000#[m]
xmax=2000#[m]
ymin=-1000#[m]
ymax=1000#[m]
dx=10#[m]
dy=10#[m]

#%% Main
x=np.arange(xmin,xmax+dx/1,dx)
y=np.arange(ymin,ymax+dy/1,dy)
X,Y=np.meshgrid(x,y)

Z=H/(1+(X/L)**2)#(Jackson-Hunt 1975 hill)


l=1/8*z0*(np.float64(L)/z0)**0.9

# #%% Output
# Output=xr.Dataset()
# Output['z']=xr.DataArray(data=Z.T,coords={'x':x,'y':y})
# Output['L']=L
# Output['H']=H
# Output.to_netcdf(os.path.join(cd,'data/Ideal_hill.nc'))

#%% Plots
fig=plt.figure(figsize=(18,3))
ax=fig.add_subplot(1,1,1)
ax.fill_between(x,x*0,Z[0,:],color='k',alpha=0.5)
ax.fill_between(x,Z[0,:],Z[0,:]+l,color='b',alpha=0.5)
ax.plot(x,Z[0,:],color='k')
