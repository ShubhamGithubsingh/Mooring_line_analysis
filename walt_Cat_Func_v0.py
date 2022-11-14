# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 23:48:31 2022

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt

def catcon_0(d,rho_w,seg_l,n_seg,spm_dry,rho_mat,d_w,X_i,g):
   
    sj_0=np.zeros(n_seg+1)
    xj_0=np.zeros(n_seg+1)
    Xj_0=np.zeros(n_seg+1)
    Zj_0=np.zeros(n_seg+1)
    lj_0=np.zeros(n_seg)
    #Initial Calculations
    #*********************************************************************************************************************
    lc=n_seg*seg_l
    X_max=np.sqrt(lc*lc-d*d)
    X0=lc-d
    xmax=X_max-10
    x=np.arange(0,xmax+0.1,0.1)
    spm_wet=spm_dry-spm_dry/rho_mat*rho_w
    
    Th=np.zeros(len(x))
    #Iteration for Th
    for i in range(len(x)):
        if x[i]==0:
            Th[i]=0
        else:
            Thl=500
            Thlc=0
            while abs(Thlc-Thl)>0.001:
                Thlc=Thl
                Th[i]=x[i]*spm_wet*g/np.arccosh(1+spm_wet*g*d/Thlc)
                Thl=Th[i]
                
    ls=np.zeros(len(x))
    X=np.zeros(len(x))
    #Th for X
    for i in range(len(x)):
        if Th[i]==0:
            ls[i]=d
            X[i]=lc-ls[i]
        else:
            ls[i]=d*np.sqrt(1+2*(Th[i]/(spm_wet*g*d)))
            X[i]=lc-ls[i]+x[i]

    for i in range(len(x)):
        X[i]=np.round(X[i],1)

    indx=np.where(X==X_i)
    indx_0=indx[0][0]
    a=Th[indx_0]/(spm_wet*g)

    x_te_init=a*np.arccosh(d/a+1)
    X_tdp=X_i-x_te_init

    count=0
    lf=0
    while lf<=X_tdp:
        sj_0[count]=lf
        Xj_0[count]=lf
        Zj_0[count]=-d
        lf=lf+seg_l
        X_seg_whole_floor=Xj_0[count]
        count=count+1

    part_sj_0_f=X_tdp-X_seg_whole_floor

    sj_0_s=0
    for i in range(count,n_seg+1):
        sj_0[i]=seg_l-part_sj_0_f+sj_0_s
        xj_0[i]=a*np.arcsinh(sj_0[i]/a)
        Zj_0[i]=-d+a*np.cosh(xj_0[i]/a)-a
        Xj_0[i]=X_tdp+xj_0[i]
        sj_0_s=sj_0_s+seg_l
    
    for i in range(0,n_seg):
        lj_0[i]=np.sqrt((Xj_0[i+1]-Xj_0[i])*(Xj_0[i+1]-Xj_0[i])+(Zj_0[i+1]-Zj_0[i])*(Zj_0[i+1]-Zj_0[i]))
    
    return(Xj_0,Zj_0,lj_0)
  
def ten_0(d,rho_w,lj_0,n_seg,spm_dry,d_w,kj,mFu,xj_c,zj_c,g,dt):
#Calculates the initial tensions in the rod links
    
    #Preliminary calculations
    #****************************************************************
    D=1.8*d_w#hydraulic dia of mooring chain
    csa_m=np.pi*D*D/4#cross sectional area of chain
    Beta_Fu=8.726547/(d-mFu*d)#parameter for sea bed upthrust
    Eta_fu=-7.25+Beta_Fu*d#parameter for sea bed upthrust
    
    #Matrix definitions
    #****************************************************************
    mj_0=np.zeros(n_seg-1)
    ejp1_0=np.zeros(n_seg-1)
    ejm1_0=np.zeros(n_seg-1)
    
    cjp1_c=np.zeros(n_seg-1)
    sjp1_c=np.zeros(n_seg-1)
    cjm1_c=np.zeros(n_seg-1)
    sjm1_c=np.zeros(n_seg-1)

    Ij_c=np.zeros(n_seg-1)
    Jj_c=np.zeros(n_seg-1)
    Kj_c=np.zeros(n_seg-1)

    Lj_c=np.zeros(n_seg-1)
    Mj_c=np.zeros(n_seg-1)
    Nj_c=np.zeros(n_seg-1)

    Pj_c=np.zeros(n_seg+1)
    Qj_c=np.zeros(n_seg+1)
    Rj_c=np.zeros(n_seg+1)
    Sj_c=np.zeros(n_seg+1)

    Wj_0=np.zeros(n_seg-1)
    Fu_c=np.zeros(n_seg-1)
  
    Xj_c=np.zeros(n_seg-1)
    Zj_c=np.zeros(n_seg-1)

    Uj_c=np.zeros(n_seg+1)
    Vj_c=np.zeros(n_seg+1)

    Ejm1_c=np.zeros(n_seg)
    Fjm1_c=np.zeros(n_seg)
    Gjm1_c=np.zeros(n_seg)
    Hjm1_c=np.zeros(n_seg)
    Psijm1_c=np.zeros(n_seg)

    alpha=np.zeros(n_seg+1)
    beta=np.zeros(n_seg+1)
    
    Ten_tent_c=np.zeros(n_seg)
    
    x_tent_n=np.zeros(n_seg+1)
    z_tent_n=np.zeros(n_seg+1)
    
    omega_tent_n=np.zeros(n_seg)

    E_tent_n=np.zeros(n_seg)
    F_tent_n=np.zeros(n_seg)
    G_tent_n=np.zeros(n_seg)

    kappa=np.zeros(n_seg+1)
    lamda=np.zeros(n_seg+1)
   
    del_Ten_c=np.zeros(n_seg)

    del_x_n=np.zeros(n_seg+1)
    del_z_n=np.zeros(n_seg+1)
    
    xj_n=np.zeros(n_seg+1)
    zj_n=np.zeros(n_seg+1)
    
    lj_n=np.zeros(n_seg)
    
    vxj_n=np.zeros(n_seg+1)
    vzj_n=np.zeros(n_seg+1)
    axj_c=np.zeros(n_seg+1)
    azj_c=np.zeros(n_seg+1)
    
    #Tentative tension at present time step
    #****************************************************************
   
    for i in range(n_seg-1):
        mj_0[i]=0.5*(spm_dry*(lj_0[i]+lj_0[i+1]))
        ejm1_0[i]=rho_w*kj*lj_0[i]*csa_m
        ejp1_0[i]=rho_w*kj*lj_0[i+1]*csa_m
        
        sjp1_c[i]=(zj_c[i+2]-zj_c[i+1])/lj_0[i+1]
        cjp1_c[i]=(xj_c[i+2]-xj_c[i+1])/lj_0[i+1]
        sjm1_c[i]=(zj_c[i+1]-zj_c[i])/lj_0[i]
        cjm1_c[i]=(xj_c[i+1]-xj_c[i])/lj_0[i]

        Ij_c[i]=mj_0[i]+0.5*(ejp1_0[i]*sjp1_c[i]*sjp1_c[i]+ejm1_0[i]*sjm1_c[i]*sjm1_c[i])
        Jj_c[i]=mj_0[i]+0.5*(ejp1_0[i]*cjp1_c[i]*cjp1_c[i]+ejm1_0[i]*cjm1_c[i]*cjm1_c[i])
        Kj_c[i]=0.5*(ejp1_0[i]*sjp1_c[i]*cjp1_c[i]+ejm1_0[i]*sjm1_c[i]*cjm1_c[i])
        
        Lj_c[i]=dt*dt*Ij_c[i]/(Ij_c[i]*Jj_c[i]-Kj_c[i]*Kj_c[i])
        Mj_c[i]=dt*dt*Jj_c[i]/(Ij_c[i]*Jj_c[i]-Kj_c[i]*Kj_c[i])
        Nj_c[i]=dt*dt*Kj_c[i]/(Ij_c[i]*Jj_c[i]-Kj_c[i]*Kj_c[i])
        
        Pj_c[i+1]=Mj_c[i]*cjm1_c[i]+Nj_c[i]*sjm1_c[i]
        Qj_c[i+1]=Nj_c[i]*cjm1_c[i]+Lj_c[i]*sjm1_c[i]
        Rj_c[i+1]=Mj_c[i]*cjp1_c[i]+Nj_c[i]*sjp1_c[i]
        Sj_c[i+1]=Nj_c[i]*cjp1_c[i]+Lj_c[i]*sjp1_c[i]
        
        Wj_0[i]=mj_0[i]*g-0.5*rho_w*g*(lj_0[i+1]*csa_m+lj_0[i]*csa_m)
    
        Fu_c[i]=Wj_0[i]*(1-np.tanh(Beta_Fu*zj_c[i+1]+Eta_fu))/2
        
        Xj_c[i]=0
        
        Zj_c[i]=-Wj_0[i]+Fu_c[i]
        
        Uj_c[i+1]=Mj_c[i]*Xj_c[i]+Nj_c[i]*Zj_c[i]
        Vj_c[i+1]=Nj_c[i]*Xj_c[i]+Lj_c[i]*Zj_c[i]
    '''
    print('mj_0','\n',mj_0,'\n')
    print('ejm1_0','\n',ejm1_0,'\n')
    print('ejp1_0','\n',ejp1_0,'\n')
    print('sjp1_c','\n',sjp1_c,'\n')
    print('cjp1_c','\n',cjp1_c,'\n')
    print('sjm1_c','\n',sjm1_c,'\n')
    print('cjm1_c','\n',cjm1_c,'\n')
    print('Ij_c','\n',Ij_c,'\n')
    print('Jj_c','\n',Jj_c,'\n')
    print('Kj_c','\n',Kj_c,'\n')
    print('Lj_c','\n',Lj_c,'\n')
    print('Mj_c','\n',Mj_c,'\n')
    print('Nj_c','\n',Nj_c,'\n')
    print('Pj_c','\n',Pj_c,'\n')
    print('Qj_c','\n',Qj_c,'\n')
    print('Rj_c','\n',Rj_c,'\n')
    print('Sj_c','\n',Sj_c,'\n')
    print('Wj_0','\n',Wj_0,'\n')
    print('Fu_c','\n',Fu_c,'\n')
    print('Xj_c','\n',Xj_c,'\n')
    print('Zj_c','\n',Zj_c,'\n')
    print('Uj_c','\n',Uj_c,'\n')
    print('Vj_c','\n',Vj_c,'\n')
    '''
    for i in range(n_seg):
        
        Ejm1_c[i]=(xj_c[i+1]-xj_c[i])*Pj_c[i]+(zj_c[i+1]-zj_c[i])*Qj_c[i]
        Fjm1_c[i]=(xj_c[i+1]-xj_c[i])*(Pj_c[i+1]+Rj_c[i])+(zj_c[i+1]-zj_c[i])*(Qj_c[i+1]+Sj_c[i])
        Gjm1_c[i]=(xj_c[i+1]-xj_c[i])*Rj_c[i+1]+(zj_c[i+1]-zj_c[i])*Sj_c[i+1]
        Hjm1_c[i]=(xj_c[i+1]-xj_c[i])*(Uj_c[i+1]-Uj_c[i])+(zj_c[i+1]-zj_c[i])*(Vj_c[i+1]-Vj_c[i])

    Ejm1_c[0]=0
    Gjm1_c[n_seg-1]=0
   
    Psijm1_c=Hjm1_c
    
    alpha[1]=1
    beta[1]=0
    for i in range(2,n_seg+1):
        if abs(Gjm1_c[i-2])>0:
            alpha[i]=(Fjm1_c[i-2]*alpha[i-1]-Ejm1_c[i-2]*alpha[i-2])/Gjm1_c[i-2]
            beta[i]=(Fjm1_c[i-2]*beta[i-1]-Ejm1_c[i-2]*beta[i-2]-Psijm1_c[i-2])/Gjm1_c[i-2]
        else:
            alpha[i]=0
            beta[i]=0
            
    den_c=(Fjm1_c[n_seg-1]*alpha[n_seg]-Ejm1_c[n_seg-1]*alpha[n_seg-1])
    if abs(den_c)<0.000001:
        T05_c=0
    else:
        T05_c=-(Fjm1_c[n_seg-1]*beta[n_seg]-Ejm1_c[n_seg-1]*beta[n_seg-1]-Hjm1_c[n_seg-1])/den_c

    Ten_tent_c[0]=T05_c
    for i in range (1,n_seg):
        Ten_tent_c[i]=alpha[i+1]*Ten_tent_c[0]+beta[i+1]
    
    
    Ten_c=Ten_tent_c
    
    print('Ejm1_c0','\n',Ejm1_c,'\n')
    print('Fjm1_c0','\n',Fjm1_c,'\n')
    print('Gjm1_c0','\n',Gjm1_c,'\n')
    print('Hjm1_c0','\n',Hjm1_c,'\n')
    '''
    print('Psijm1_c0','\n',Psijm1_c,'\n')
    print('aplha_0','\n',alpha,'\n')
    print('beta_0','\n',beta,'\n')
    print('den_c','\n',den_c,'\n')
    print('T05_c0','\n',T05_c,'\n')
    print('Ten_tent_c0','\n',Ten_tent_c,'\n')
    '''
    omega_tent_n=np.ones(n_seg)
    
    itr_count=0
    
    kappa[1]=1
    lamda[1]=0
    
    while np.sum(abs(omega_tent_n))>0.01*n_seg:
        if itr_count>1000:
            print('1000 iterations exceeded!')
            exit()
    
        for i in range(n_seg+1):
            if i==0:
                x_tent_n[i]=xj_c[i]
                z_tent_n[i]=zj_c[i]
            elif i==n_seg:
                x_tent_n[i]=xj_c[i]+Uj_c[i]
                z_tent_n[i]=zj_c[i]+Vj_c[i]
            else:
                x_tent_n[i]=xj_c[i]+0.5*(-Pj_c[i]*Ten_c[i-1]+Rj_c[i]*Ten_c[i]+Uj_c[i])
                z_tent_n[i]=zj_c[i]+0.5*(-Qj_c[i]*Ten_c[i-1]+Sj_c[i]*Ten_c[i]+Vj_c[i])
                
        for i in range(n_seg):
            omega_tent_n[i]=0.5*((x_tent_n[i+1]-x_tent_n[i])*(x_tent_n[i+1]-x_tent_n[i])
                                +(z_tent_n[i+1]-z_tent_n[i])*(z_tent_n[i+1]-z_tent_n[i])-lj_0[i]*lj_0[i])
        
            E_tent_n[i]=(x_tent_n[i+1]-x_tent_n[i])*Pj_c[i]+(z_tent_n[i+1]-z_tent_n[i])*Qj_c[i]
            F_tent_n[i]=(x_tent_n[i+1]-x_tent_n[i])*(Pj_c[i+1]+Rj_c[i])+(z_tent_n[i+1]-z_tent_n[i])*(Qj_c[i+1]+Sj_c[i])
            G_tent_n[i]=(x_tent_n[i+1]-x_tent_n[i])*Rj_c[i+1]+(z_tent_n[i+1]-z_tent_n[i])*Sj_c[i+1]
        
        for i in range(2,n_seg+1):
            if abs(G_tent_n[i-2])>0:
                kappa[i]=(F_tent_n[i-2]*kappa[i-1]-E_tent_n[i-2]*kappa[i-2])/G_tent_n[i-2]
                lamda[i]=(F_tent_n[i-2]*lamda[i-1]-E_tent_n[i-2]*lamda[i-2]-omega_tent_n[i-2])/G_tent_n[i-2]
            else:
                kappa[i]=0
                lamda[i]=0

        den_n=(F_tent_n[n_seg-1]*kappa[n_seg]-E_tent_n[n_seg-1]*kappa[n_seg-1])
        if abs(den_n)<0.000001:
            del_T05_c=0
        else:
            del_T05_c=-(F_tent_n[n_seg-1]*lamda[n_seg]-E_tent_n[n_seg-1]*lamda[n_seg-1]-omega_tent_n[n_seg-1])/den_n

        del_Ten_c[0]=del_T05_c
        for i in range (1,n_seg):
            del_Ten_c[i]=kappa[i+1]*del_Ten_c[0]+lamda[i+1]
        
        Ten_c=Ten_c+del_Ten_c
              
        itr_count=itr_count+1
    
    for i in range(1,n_seg):
        del_x_n[i]=0.5*(-Pj_c[i]*del_Ten_c[i-1]+Rj_c[i]*del_Ten_c[i])
        del_z_n[i]=0.5*(-Qj_c[i]*del_Ten_c[i-1]+Sj_c[i]*del_Ten_c[i])

    xj_n=x_tent_n+del_x_n
    zj_n=z_tent_n+del_z_n
           
    for i in range(n_seg):
        lj_n[i]=np.sqrt((xj_n[i+1]-xj_n[i])*(xj_n[i+1]-xj_n[i])+(zj_n[i+1]-zj_n[i])*(zj_n[i+1]-zj_n[i]))
        
    for i in range(n_seg+1):
        vxj_n[i]=(xj_n[i]-xj_c[i])/dt
        vzj_n[i]=(zj_n[i]-zj_c[i])/dt
        axj_c[i]=(xj_n[i]-xj_c[i])/(dt*dt)
        azj_c[i]=(zj_n[i]-zj_c[i])/(dt*dt)
    ''' 
    print('x_tent_n0','\n',x_tent_n,'\n')
    print('z_tent_n0','\n',z_tent_n,'\n')
    print('omega_tent_n0','\n',omega_tent_n,'\n')
    print('E_tent_n','\n',E_tent_n,'\n')
    print('F_tent_n','\n',F_tent_n,'\n')
    print('G_tent_n','\n',G_tent_n,'\n')
    print('kappa_0','\n',kappa,'\n')
    print('lamda_0','\n',lamda,'\n')
    print('den_n','\n',den_c,'\n')
    print('del_T05_c0','\n',T05_c,'\n')
    print('del_Ten_c0','\n',del_Ten_c,'\n')
    '''
    print('Ten_c0','\n',Ten_c,'\n')
    '''
    print('del_x_n0','\n',del_x_n,'\n')
    print('del_z_n0','\n',del_z_n,'\n')
    print('itr_count0','\n',itr_count,'\n')
    '''
    print('xj_n0','\n',xj_n,'\n')
    print('zj_n0','\n',zj_n,'\n')
    print('vxj_n0','\n',vxj_n,'\n')
    print('vzj_n0','\n',vzj_n,'\n')
    print('axj_c0','\n',axj_c,'\n')
    print('azj_c0','\n',azj_c,'\n')
    '''
    print('lj_n0','\n',lj_n,'\n')
    '''
    return(xj_c,zj_c,Ten_c,Ejm1_c,Fjm1_c,Gjm1_c,Hjm1_c,xj_n,zj_n,lj_n,vxj_n,vzj_n,axj_c,azj_c)

def ten(d,rho_w,c,lj_0,n_seg,spm_dry,d_w,kj,cd,mFu,xj_c,zj_c,vxj_c,vzj_c,axj_p,azj_p,xj_p,zj_p,Ten_p,Ejm1_p,Fjm1_p,Gjm1_p,Hjm1_p,g,dt,U_x):
#Calculates the initial tensions in the rod links
    
    #Preliminary calculations
    #****************************************************************
    D=1.8*d_w#hydraulic dia of mooring chain
    csa_m=np.pi*D*D/4#cross sectional area of chain
    Beta_Fu=8.726547/(d-mFu*d)#parameter for sea bed upthrust
    Eta_fu=-7.25+Beta_Fu*d#parameter for sea bed upthrust
    
    #Matrix definitions
    #****************************************************************
    mj_0=np.zeros(n_seg-1)
    ejp1_0=np.zeros(n_seg-1)
    ejm1_0=np.zeros(n_seg-1)
    
    cjp1_c=np.zeros(n_seg-1)
    sjp1_c=np.zeros(n_seg-1)
    cjm1_c=np.zeros(n_seg-1)
    sjm1_c=np.zeros(n_seg-1)

    Ij_c=np.zeros(n_seg-1)
    Jj_c=np.zeros(n_seg-1)
    Kj_c=np.zeros(n_seg-1)

    Lj_c=np.zeros(n_seg-1)
    Mj_c=np.zeros(n_seg-1)
    Nj_c=np.zeros(n_seg-1)

    Pj_c=np.zeros(n_seg+1)
    Qj_c=np.zeros(n_seg+1)
    Rj_c=np.zeros(n_seg+1)
    Sj_c=np.zeros(n_seg+1)

    Wj_0=np.zeros(n_seg-1)
    Fu_c=np.zeros(n_seg-1)
    
    qjp1_c=np.zeros(n_seg-1)
    qjm1_c=np.zeros(n_seg-1)

    fjp1_c=np.zeros(n_seg-1)
    fjm1_c=np.zeros(n_seg-1)

    Djp1_c=np.zeros(n_seg-1)
    Djm1_c=np.zeros(n_seg-1)
  
    Xj_c=np.zeros(n_seg-1)
    Zj_c=np.zeros(n_seg-1)

    Uj_c=np.zeros(n_seg+1)
    Vj_c=np.zeros(n_seg+1)

    Ejm1_c=np.zeros(n_seg)
    Fjm1_c=np.zeros(n_seg)
    Gjm1_c=np.zeros(n_seg)
    Hjm1_c=np.zeros(n_seg)
    Psijm1_c=np.zeros(n_seg)

    alpha=np.zeros(n_seg+1)
    beta=np.zeros(n_seg+1)
    
    Ten_tent_c=np.zeros(n_seg)
    
    x_tent_n=np.zeros(n_seg+1)
    z_tent_n=np.zeros(n_seg+1)
    
    omega_tent_n=np.zeros(n_seg)

    E_tent_n=np.zeros(n_seg)
    F_tent_n=np.zeros(n_seg)
    G_tent_n=np.zeros(n_seg)

    kappa=np.zeros(n_seg+1)
    lamda=np.zeros(n_seg+1)
   
    del_Ten_c=np.zeros(n_seg)

    del_x_n=np.zeros(n_seg+1)
    del_z_n=np.zeros(n_seg+1)
    
    xj_n=np.zeros(n_seg+1)
    zj_n=np.zeros(n_seg+1)
    
    lj_n=np.zeros(n_seg)
    
    vxj_n=np.zeros(n_seg+1)
    vzj_n=np.zeros(n_seg+1)
    axj_c=np.zeros(n_seg+1)
    azj_c=np.zeros(n_seg+1)
    
    #Tentative tension at present time step
    #****************************************************************
   
    for i in range(n_seg-1):
        mj_0[i]=0.5*(spm_dry*(lj_0[i]+lj_0[i+1]))
        ejm1_0[i]=rho_w*kj*lj_0[i]*csa_m
        ejp1_0[i]=rho_w*kj*lj_0[i+1]*csa_m
        
        sjp1_c[i]=(zj_c[i+2]-zj_c[i+1])/lj_0[i+1]
        cjp1_c[i]=(xj_c[i+2]-xj_c[i+1])/lj_0[i+1]
        sjm1_c[i]=(zj_c[i+1]-zj_c[i])/lj_0[i]
        cjm1_c[i]=(xj_c[i+1]-xj_c[i])/lj_0[i]

        Ij_c[i]=mj_0[i]+0.5*(ejp1_0[i]*sjp1_c[i]*sjp1_c[i]+ejm1_0[i]*sjm1_c[i]*sjm1_c[i])
        Jj_c[i]=mj_0[i]+0.5*(ejp1_0[i]*cjp1_c[i]*cjp1_c[i]+ejm1_0[i]*cjm1_c[i]*cjm1_c[i])
        Kj_c[i]=0.5*(ejp1_0[i]*sjp1_c[i]*cjp1_c[i]+ejm1_0[i]*sjm1_c[i]*cjm1_c[i])
        
        Lj_c[i]=dt*dt*Ij_c[i]/(Ij_c[i]*Jj_c[i]-Kj_c[i]*Kj_c[i])
        Mj_c[i]=dt*dt*Jj_c[i]/(Ij_c[i]*Jj_c[i]-Kj_c[i]*Kj_c[i])
        Nj_c[i]=dt*dt*Kj_c[i]/(Ij_c[i]*Jj_c[i]-Kj_c[i]*Kj_c[i])
        
        Pj_c[i+1]=Mj_c[i]*cjm1_c[i]+Nj_c[i]*sjm1_c[i]
        Qj_c[i+1]=Nj_c[i]*cjm1_c[i]+Lj_c[i]*sjm1_c[i]
        Rj_c[i+1]=Mj_c[i]*cjp1_c[i]+Nj_c[i]*sjp1_c[i]
        Sj_c[i+1]=Nj_c[i]*cjp1_c[i]+Lj_c[i]*sjp1_c[i]
        
        qjp1_c[i]=-0.5*((vxj_c[i+2]-c)+(vxj_c[i+1]-c))*sjp1_c[i]+0.5*(vzj_c[i+2]+vzj_c[i+1])*cjp1_c[i]
        qjm1_c[i]=-0.5*((vxj_c[i]-c)+(vxj_c[i+1]-c))*sjm1_c[i]+0.5*(vzj_c[i]+vzj_c[i+1])*cjm1_c[i]
        
        fjp1_c[i]=0.5*rho_w*cd*lj_0[i+1]*D
        fjm1_c[i]=0.5*rho_w*cd*lj_0[i]*D
        
        Djp1_c[i]=-fjp1_c[i]*qjp1_c[i]*abs(qjp1_c[i])
        Djm1_c[i]=-fjm1_c[i]*qjm1_c[i]*abs(qjm1_c[i])
        
        Wj_0[i]=mj_0[i]*g-0.5*rho_w*g*(lj_0[i+1]*csa_m+lj_0[i]*csa_m)
    
        Fu_c[i]=Wj_0[i]*(1-np.tanh(Beta_Fu*zj_c[i+1]+Eta_fu))/2
        
        Xj_c[i]=-0.5*(Djp1_c[i]*sjp1_c[i]+Djm1_c[i]*sjm1_c[i])
        
        Zj_c[i]=-Wj_0[i]+Fu_c[i]
        
        Uj_c[i+1]=Mj_c[i]*Xj_c[i]+Nj_c[i]*Zj_c[i]
        Vj_c[i+1]=Nj_c[i]*Xj_c[i]+Lj_c[i]*Zj_c[i]
    Uj_c[n_seg]=U_x
    '''
    print('mj_0','\n',mj_0,'\n')
    print('ejm1_0','\n',ejm1_0,'\n')
    print('ejp1_0','\n',ejp1_0,'\n')
    print('sjp1_c','\n',sjp1_c,'\n')
    print('cjp1_c','\n',cjp1_c,'\n')
    print('sjm1_c','\n',sjm1_c,'\n')
    print('cjm1_c','\n',cjm1_c,'\n')
    print('Ij_c','\n',Ij_c,'\n')
    print('Jj_c','\n',Jj_c,'\n')
    print('Kj_c','\n',Kj_c,'\n')
    print('Lj_c','\n',Lj_c,'\n')
    print('Mj_c','\n',Mj_c,'\n')
    print('Nj_c','\n',Nj_c,'\n')
    print('Pj_c','\n',Pj_c,'\n')
    print('Qj_c','\n',Qj_c,'\n')
    print('Rj_c','\n',Rj_c,'\n')
    print('Sj_c','\n',Sj_c,'\n')
    print('Djm1_c','\n',Djm1_c,'\n')
    print('Djp1_c','\n',Djp1_c,'\n')
    print('Wj_0','\n',Wj_0,'\n')
    print('Fu_c','\n',Fu_c,'\n')
    print('Xj_c','\n',Xj_c,'\n')
    print('Zj_c','\n',Zj_c,'\n')
    print('Uj_c','\n',Uj_c,'\n')
    print('Vj_c','\n',Vj_c,'\n')
    '''
    for i in range(n_seg):
        
        Ejm1_c[i]=(xj_c[i+1]-xj_c[i])*Pj_c[i]+(zj_c[i+1]-zj_c[i])*Qj_c[i]
        Fjm1_c[i]=(xj_c[i+1]-xj_c[i])*(Pj_c[i+1]+Rj_c[i])+(zj_c[i+1]-zj_c[i])*(Qj_c[i+1]+Sj_c[i])
        Gjm1_c[i]=(xj_c[i+1]-xj_c[i])*Rj_c[i+1]+(zj_c[i+1]-zj_c[i])*Sj_c[i+1]
        Hjm1_c[i]=(xj_c[i+1]-xj_c[i])*(Uj_c[i+1]-Uj_c[i])+(zj_c[i+1]-zj_c[i])*(Vj_c[i+1]-Vj_c[i])

    Ejm1_c[0]=0
    Gjm1_c[n_seg-1]=0
   
    for i in range(n_seg):
        if i==0:
            Psijm1_c[i]=(-Fjm1_p[i]*Ten_p[i]+Gjm1_p[i]*Ten_p[i+1]+Hjm1_p[i]+Hjm1_c[i]
                    +2*((xj_c[i+1]-xj_p[i+1])-(xj_c[i]-xj_p[i]))*((xj_c[i+1]-xj_p[i+1])-(xj_c[i]-xj_p[i]))
                    +2*((zj_c[i+1]-zj_p[i+1])-(zj_c[i]-zj_p[i]))*((zj_c[i+1]-zj_p[i+1])-(zj_c[i]-zj_p[i])))
        if i==n_seg-1:
            Psijm1_c[i]=(Ejm1_p[i]*Ten_p[i-1]-Fjm1_p[i]*Ten_p[i]+Hjm1_p[i]+Hjm1_c[i]
                    +2*((xj_c[i+1]-xj_p[i+1])-(xj_c[i]-xj_p[i]))*((xj_c[i+1]-xj_p[i+1])-(xj_c[i]-xj_p[i]))
                    +2*((zj_c[i+1]-zj_p[i+1])-(zj_c[i]-zj_p[i]))*((zj_c[i+1]-zj_p[i+1])-(zj_c[i]-zj_p[i])))
        else:   
            Psijm1_c[i]=(Ejm1_p[i]*Ten_p[i-1]-Fjm1_p[i]*Ten_p[i]+Gjm1_p[i]*Ten_p[i+1]+Hjm1_p[i]+Hjm1_c[i]
                    +2*((xj_c[i+1]-xj_p[i+1])-(xj_c[i]-xj_p[i]))*((xj_c[i+1]-xj_p[i+1])-(xj_c[i]-xj_p[i]))
                    +2*((zj_c[i+1]-zj_p[i+1])-(zj_c[i]-zj_p[i]))*((zj_c[i+1]-zj_p[i+1])-(zj_c[i]-zj_p[i])))
    
    alpha[1]=1
    beta[1]=0
    for i in range(2,n_seg+1):
        if abs(Gjm1_c[i-2])>0:
            alpha[i]=(Fjm1_c[i-2]*alpha[i-1]-Ejm1_c[i-2]*alpha[i-2])/Gjm1_c[i-2]
            beta[i]=(Fjm1_c[i-2]*beta[i-1]-Ejm1_c[i-2]*beta[i-2]-Psijm1_c[i-2])/Gjm1_c[i-2]
        else:
            alpha[i]=0
            beta[i]=0
            
    den_c=(Fjm1_c[n_seg-1]*alpha[n_seg]-Ejm1_c[n_seg-1]*alpha[n_seg-1])
    if abs(den_c)<0.000001:
        T05_c=0
    else:
        T05_c=-(Fjm1_c[n_seg-1]*beta[n_seg]-Ejm1_c[n_seg-1]*beta[n_seg-1]-Hjm1_c[n_seg-1])/den_c

    Ten_tent_c[0]=T05_c
    for i in range (1,n_seg):
        Ten_tent_c[i]=alpha[i+1]*Ten_tent_c[0]+beta[i+1]
    
    
    Ten_c=Ten_tent_c
    '''
    print('Ejm1_c','\n',Ejm1_c,'\n')
    print('Fjm1_c','\n',Fjm1_c,'\n')
    print('Gjm1_c','\n',Gjm1_c,'\n')
    print('Hjm1_c','\n',Hjm1_c,'\n')
    print('Psijm1_c','\n',Psijm1_c,'\n')
    print('aplha','\n',alpha,'\n')
    print('beta','\n',beta,'\n')
    print('den_c','\n',den_c,'\n')
    print('T05_c','\n',T05_c,'\n')
    print('Ten_tent_c','\n',Ten_tent_c,'\n')
    '''
    omega_tent_n=np.ones(n_seg)
    
    itr_count=0
    
    kappa[1]=1
    lamda[1]=0
    
    while np.sum(abs(omega_tent_n))>0.01*n_seg:
        if itr_count>1000:
            print('1000 iterations exceeded!')
            exit()
    
        for i in range(n_seg+1):
            if i==0:
                x_tent_n[i]=xj_c[i]
                z_tent_n[i]=zj_c[i]
            elif i==n_seg:
                x_tent_n[i]=xj_c[i]+Uj_c[i]
                z_tent_n[i]=zj_c[i]+Vj_c[i]
            else:
                x_tent_n[i]=2*xj_c[i]-xj_p[i] -Pj_c[i] * Ten_c[i - 1] + Rj_c[i] * Ten_c[i] + Uj_c[i]
                z_tent_n[i]=2*zj_c[i]-zj_p[i] -Qj_c[i] * Ten_c[i - 1] + Sj_c[i] * Ten_c[i] + Vj_c[i]
                
        for i in range(n_seg):
            omega_tent_n[i]=0.5*((x_tent_n[i+1]-x_tent_n[i])*(x_tent_n[i+1]-x_tent_n[i])
                                +(z_tent_n[i+1]-z_tent_n[i])*(z_tent_n[i+1]-z_tent_n[i])-lj_0[i]*lj_0[i])
        
            E_tent_n[i]=(x_tent_n[i+1]-x_tent_n[i])*Pj_c[i]+(z_tent_n[i+1]-z_tent_n[i])*Qj_c[i]
            F_tent_n[i]=(x_tent_n[i+1]-x_tent_n[i])*(Pj_c[i+1]+Rj_c[i])+(z_tent_n[i+1]-z_tent_n[i])*(Qj_c[i+1]+Sj_c[i])
            G_tent_n[i]=(x_tent_n[i+1]-x_tent_n[i])*Rj_c[i+1]+(z_tent_n[i+1]-z_tent_n[i])*Sj_c[i+1]
        
        for i in range(2,n_seg+1):
            if abs(G_tent_n[i-2])>0:
                kappa[i]=(F_tent_n[i-2]*kappa[i-1]-E_tent_n[i-2]*kappa[i-2])/G_tent_n[i-2]
                lamda[i]=(F_tent_n[i-2]*lamda[i-1]-E_tent_n[i-2]*lamda[i-2]-omega_tent_n[i-2])/G_tent_n[i-2]
            else:
                kappa[i]=0
                lamda[i]=0

        den_n=(F_tent_n[n_seg-1]*kappa[n_seg]-E_tent_n[n_seg-1]*kappa[n_seg-1])
        if abs(den_n)<0.000001:
            del_T05_c=0
        else:
            del_T05_c=-(F_tent_n[n_seg-1]*lamda[n_seg]-E_tent_n[n_seg-1]*lamda[n_seg-1]-omega_tent_n[n_seg-1])/den_n

        del_Ten_c[0]=del_T05_c
        for i in range (1,n_seg):
            del_Ten_c[i]=kappa[i+1]*del_Ten_c[0]+lamda[i+1]
        
        Ten_c=Ten_c+del_Ten_c
               
        itr_count=itr_count+1
    
    for i in range(1,n_seg):
        del_x_n[i]=-Pj_c[i]*del_Ten_c[i-1]+Rj_c[i]*del_Ten_c[i]
        del_z_n[i]=-Qj_c[i]*del_Ten_c[i-1]+Sj_c[i]*del_Ten_c[i]

    xj_n=x_tent_n+del_x_n
    zj_n=z_tent_n+del_z_n
           
    for i in range(n_seg):
        lj_n[i]=np.sqrt((xj_n[i+1]-xj_n[i])*(xj_n[i+1]-xj_n[i])+(zj_n[i+1]-zj_n[i])*(zj_n[i+1]-zj_n[i]))
        
    for i in range(n_seg+1):
        vxj_n[i]=(xj_n[i]-xj_c[i])/dt
        vzj_n[i]=(zj_n[i]-zj_c[i])/dt
        axj_c[i]=(xj_n[i]-2*xj_c[i]+xj_p[i])/(dt*dt)
        azj_c[i]=(zj_n[i]-2*zj_c[i]+zj_p[i])/(dt*dt)
    '''
    print('x_tent_n','\n',x_tent_n,'\n')
    print('z_tent_n','\n',z_tent_n,'\n')
    print('omega_tent_n','\n',omega_tent_n,'\n')
    print('E_tent_n','\n',E_tent_n,'\n')
    print('F_tent_n','\n',F_tent_n,'\n')
    print('G_tent_n','\n',G_tent_n,'\n')
    print('kappa','\n',kappa,'\n')
    print('lamda','\n',lamda,'\n')
    print('den_n','\n',den_c,'\n')
    print('del_T05_c','\n',T05_c,'\n')
    print('del_Ten_c','\n',del_Ten_c,'\n')
    print('Ten_c','\n',Ten_c,'\n')
    print('del_x_n','\n',del_x_n,'\n')
    print('del_z_n','\n',del_z_n,'\n')
    print('itr_count','\n',itr_count,'\n')
    print('xj_n','\n',xj_n,'\n')
    print('zj_n','\n',zj_n,'\n')
    print('lj_n','\n',lj_n,'\n')
    '''
    return(xj_c,zj_c,Ten_c,Ejm1_c,Fjm1_c,Gjm1_c,Hjm1_c,xj_n,zj_n,lj_n,vxj_n,vzj_n,axj_c,azj_c)
