# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 23:47:25 2022

@author: user
"""
import numpy as np
import matplotlib.pyplot as plt
import walt_Cat_Func_v0

#Parameters
#*********************************************************************************************************************
#Environment
rho_w=1025#density of water
d=50 #water depth
x_cur_max=0 #uniform current in x direction
#Mooring
rho_mat=7800 #density of chain material
seg_l=10#mooring line length
n_seg=int(10) #number of links
spm_dry=10#linear density of chain kg/m
d_w=0.022#wire dia of mooring chain
X_i=55#initial position of top end of mooring
cd=1#coefficient of drag
kj=1#added mass coefficient
mp=0#mass of point load (none considered)
vp=0#volume of point load (non considered)
#Miscellaneous
sim_time=60#simulation time in seconds
ramp_time=10#ramping of current to full value
mov_time=60#time for which top end moves
dt=0.1#time interval
g=9.81#acceleration of gravity
mFu=0.999999#factor of depth at which  vertical force becomes 5% of nodal weight
T_s=np.arange(0,sim_time+dt,dt)
T_hist=np.zeros(len(T_s))
print(len(T_s))
#Initial geometric configuration
#*********************************************************************************************************************
(xj_0,zj_0,lj_0)=walt_Cat_Func_v0.catcon_0(d,rho_w,seg_l,n_seg,spm_dry,rho_mat,d_w,X_i,g)
print('xj_0','\n',xj_0,'\n','zj_0','\n',zj_0,'\n','lj_0','\n',lj_0,'\n')

#Initial tensions
#*********************************************************************************************************************

(xj_c,zj_c,Ten_c,Ejm1_c,Fjm1_c,Gjm1_c,Hjm1_c,xj_n,zj_n,lj_n,vxj_n,vzj_n,axj_c,azj_c)=walt_Cat_Func_v0.ten_0(d,rho_w,lj_0,
                                        n_seg,spm_dry,d_w,kj,mFu,xj_0,zj_0,g,dt)
i=1
T_hist[i]=Ten_c[n_seg-1]
plt.figure(figsize=(6,5))
plt.figure(0)
plt.plot(xj_n,zj_n)
'''
print('xj_c0','\n',xj_c,'\n')
print('zj_c0','\n',zj_c,'\n')
print('Ten_c0','\n',Ten_c,'\n')
print('Ejm1_c0','\n',Ejm1_c,'\n')
print('Fjm1_c0','\n',Fjm1_c,'\n')
print('Gjm1_c0','\n',Gjm1_c,'\n')
print('Hjm1_c0','\n',Hjm1_c,'\n')
print('xj_n0','\n',xj_n,'\n')
print('zj_n0','\n',zj_n,'\n')
print('lj_n0','\n',lj_n,'\n')
print('vxj_n0','\n',vxj_n,'\n')
print('vzj_n0','\n',vzj_n,'\n')
print('axj_c0','\n',axj_c,'\n')
print('azj_c0','\n',azj_c,'\n')
'''
#Time Stepping
#*********************************************************************************************************************
time=dt
i=2
while time<=sim_time:
    c=x_cur_max
    if time<ramp_time:
        c=c*time/ramp_time
    if time<=mov_time:
        U_x=dt*dt*4
    else:
        U_x=0
    xj_p=xj_c
    zj_p=zj_c
    xj_c=xj_n
    zj_c=zj_n
    Ten_p=Ten_c
    Ejm1_p=Ejm1_c
    Fjm1_p=Fjm1_c
    Gjm1_p=Gjm1_c
    Hjm1_p=Hjm1_c
    vxj_c=vxj_n
    vzj_c=vzj_n
    axj_p=axj_c
    azj_p=azj_c
    (xj_c,zj_c,Ten_c,Ejm1_c,Fjm1_c,Gjm1_c,Hjm1_c,xj_n,zj_n,lj_n,vxj_n,vzj_n,axj_c,azj_c)=walt_Cat_Func_v0.ten(d,rho_w,c,
    lj_0,n_seg,spm_dry,d_w,kj,cd,mFu,xj_c,zj_c,vxj_c,vzj_c,axj_p,azj_p,xj_p,zj_p,Ten_p,Ejm1_p,Fjm1_p,Gjm1_p,Hjm1_p,g,dt,U_x)
    T_hist[i]=Ten_c[n_seg-1]
    time=time+dt
    i=i+1
'''
print('xj_n','\n',xj_n,'\n')
print('zj_n','\n',zj_n,'\n')
print('Ten_c','\n',Ten_c,'\n')
print('Ejm1_c','\n',Ejm1_c,'\n')
print('Fjm1_c','\n',Fjm1_c,'\n')
print('Gjm1_c','\n',Gjm1_c,'\n')
print('Hjm1_c','\n',Hjm1_c,'\n')
print('lj_n','\n',lj_n,'\n')
print('vxj_n','\n',vxj_n,'\n')
print('vzj_n','\n',vzj_n,'\n')
print('axj_c','\n',axj_c,'\n')
print('azj_c','\n',azj_c,'\n')
'''
print(Ten_c)
print(xj_c)
plt.figure(0)
plt.plot(xj_n,zj_n)
plt.figure(1)
plt.plot(T_s, T_hist)
plt.show()