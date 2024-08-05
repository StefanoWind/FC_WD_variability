# -*- coding: utf-8 -*-
"""
Functions for flow in complex terrain
"""
import numpy as np

def speedup(x,z,tapering,sigma,z0):
    reals=~np.isnan(z)
    x_sel=x[reals]
    z_sel=z[reals]
    lf=np.polyfit(x_sel, z_sel,1)
    z_det=z_sel-lf[0]*x_sel-lf[1]
    
    #full tapering matrix
    if not hasattr(speedup, "T"):
        T=np.eye(len(x))
        for i in range(len(x)):
            if x[i]-x[0]<tapering:
                T[i,i]=3*((x[i]-x[0])/tapering)**2-2*((x[i]-x[0])/tapering)**3
            if x[-1]-x[i]<tapering:
                T[i,i]=3*((x[-1]-x[i])/tapering)**2-2*((x[-1]-x[i])/tapering)**3
        speedup.T=T
    
    #full smoothing matrix
    if not hasattr(speedup, "M"):
        M=np.zeros((len(x),len(x)))
        for j in range(len(x)):
            M[j,:]=np.exp(-(x-x[j])**2/(2*sigma**2))
            M[j,:]=M[j,:]/np.sum(M[j,:])
        speedup.M=M
    
    #tapering matrix (if nans are present)
    if np.sum(reals)<len(x):
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
    else:
        T=speedup.T
        M=speedup.M

    z_smooth=np.matmul(M,z_det)#smoothed topography
    z_taper=np.matmul(T,z_smooth-np.nanmin(z_smooth))#tapering
    H=np.nanmax(z_taper)#height
    f=z_taper/H#non-dimensional hill function
    
    #lengthscales
    i1=np.where(f>0.5)[0][0]
    i2=np.where(f>0.5)[0][-1]
    L=(x_sel[i2]-x_sel[i1])/2
    l=1/8*z0*(L/z0)**0.9
    
    #weighting function
    dx=x_sel[1]-x_sel[0]
    X2,X3=np.meshgrid(x_sel,x_sel)
    X_diff=(X2-X3)
    X_diff[np.abs(X_diff)<(dx+0.0)/10]=10**10
    
    df_dx=np.concatenate([[(f[1]-f[0])/dx],(f[2:]-f[:-2])/(2*dx), [(f[-1]-f[-2])/dx]])
    dF_dXX=np.tile(df_dx,(len(x_sel),1)).T
    I=dF_dXX*L/X_diff
    I[np.isnan(I)]=0
    
    DS=np.zeros(len(x))+np.nan
    DS[reals]=H/L*np.log(L/z0)**2/np.log(l/z0)**2*np.sum(I,axis=0)/np.pi*dx*100
    
    f_all=np.zeros(len(x))+np.nan
    f_all[reals]=f
    return DS,f_all,H,L,l