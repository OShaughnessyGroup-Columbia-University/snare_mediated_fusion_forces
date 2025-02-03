# IMPORT STUFF
###########################################################
from __future__ import division
import itertools
from random import randint
import math
#import os
#import glob
#import re
#import pandas as pd
import numpy as np
#from lxml import etree
#from io import StringIO
#from lipid_classes import *
#from potentials import *
#from simple_bilayer import create_bilayer1
import sys

import hoomd
import hoomd.md

#import hoomd.deprecated
hoomd.context.initialize()

##### Define rigid bodies in the system

##### Define vesicle #####
def vesicle(radius,z_center,rho,beads_per_lipid,r_tail,r_head,n_tmd,tmd_theta):
    hmem = 2 * ((beads_per_lipid - 1) * 2 * r_tail + 2 * r_head) - 2*r_head
    dpt=np.sqrt(1/rho) # distance between lipids along phi and theta direction
    z_c = z_center
    pos_x = []
    pos_y = []
    pos_z = []
    tmd_x=[]
    tmd_y=[]
    tmd_z=[]
    tmd_orien=[]
    r1=radius+r_head # inner leaflet head bead shell radius
    r2=r1+hmem # outer leaflet head bead shell radius
    n_theta1=int(np.pi*r1/dpt)+1
    n_theta2=int(np.pi*r2/dpt)+1
    theta_head1=np.linspace(0,np.pi,n_theta1)
    theta_head1=((theta_head1+np.roll(theta_head1,-1))/2)[0:-1]
    theta_head2=np.linspace(0,np.pi,n_theta2)
    theta_head2=((theta_head2+np.roll(theta_head2,-1))/2)[0:-1]
    xh1,yh1,zh1,xtmd_in,ytmd_in,ztmd_in=place_shell(r1,theta_head1,dpt,n_tmd,tmd_theta,False)
    xh2,yh2,zh2,xtmd_out,ytmd_out,ztmd_out=place_shell(r2,theta_head2,dpt,n_tmd,tmd_theta,True)
    xh=np.concatenate((xh1,xh2))
    yh=np.concatenate((yh1,yh2))
    zh=np.concatenate((zh1,zh2))
    xtmd=np.concatenate((xtmd_in,xtmd_out))
    ytmd=np.concatenate((ytmd_in,ytmd_out))
    ztmd=np.concatenate((ztmd_in,ztmd_out))
    if n_tmd >0:
        tmd_d_phi=2*np.pi/n_tmd
    for (i,x) in enumerate(xh):
        yhh=yh[i]
        zhh=zh[i]
        r=np.sqrt((x**2)+(yhh**2)+(zhh**2))
        if r < radius+0.5:
            nrm=[x,yhh,zhh]
        else:
            nrm=[-x,-yhh,-zhh]
        x,y,z=place_tail(x,yhh,zhh,nrm,beads_per_lipid)
        z=[zl+z_c for zl in z]
        pos_x.extend(x)
        pos_y.extend(y)
        pos_z.extend(z)
    for i in range(n_tmd):
        x=xtmd[i]
        y=ytmd[i]
        z=ztmd[i]
        nrm=np.array([x,y,z])
        nrm=nrm/(np.sqrt((nrm[0]**2)+(nrm[1]**2)+(nrm[2]**2)))
        tmd_orien_i = np.array([1,0,0])
        rot_axis = np.cross(nrm,tmd_orien_i)
        rot_axis=rot_axis/(np.sqrt((rot_axis[0]**2)+(rot_axis[1]**2)+(rot_axis[2]**2)))
        nrm_dot_tmd_orien = np.dot(nrm,tmd_orien_i)
        rot_ang = -np.arccos(nrm_dot_tmd_orien)
        tmd_x.append(x)
        tmd_y.append(y)
        tmd_z.append(z+z_c)
        ang=tmd_theta
        orientation=[np.cos(rot_ang/2),rot_axis[0]*np.sin(rot_ang/2),rot_axis[1]*np.sin(rot_ang/2),rot_axis[2]*np.sin(rot_ang/2)]#[np.cos(ang/2),0,0,np.sin(ang/2)]#[0,x/r,y/r,z/r]
        tmd_orien.append(orientation)
    print('num beads %.1f %.1f %.1f' %(len(pos_x),len(pos_y),len(pos_z)))
    return(np.array(pos_x),np.array(pos_y),np.array(pos_z),np.array(tmd_x),np.array(tmd_y),np.array(tmd_z),tmd_orien)

def place_tail(x,y,z,nrm,beads_per_lipid):
    nrm=nrm/(np.sqrt((nrm[0]**2)+(nrm[1]**2)+(nrm[2]**2)))
    xl=[]
    yl=[]
    zl=[]
    xl.append(x)
    yl.append(y)
    zl.append(z)
    for i in range(beads_per_lipid-1):
        xl.append(xl[i]+nrm[0])
        yl.append(yl[i]+nrm[1])
        zl.append(zl[i]+nrm[2])
    return(xl,yl,zl)

def place_shell(r,theta,dpt,n_tmd,tmd_theta,ec_lef):
    theta_head1=theta
    r1=r
    xh=[]
    yh=[]
    zh=[]
    xtmd=[]
    ytmd=[]
    ztmd=[]
    tmd_phi=[]
    d=3
    z_tmd = (r1-d)*math.cos(tmd_theta)
    r_tmd = (r1-d)*np.sin(tmd_theta)
    dt_lip = theta_head1[1]-theta_head1[0]
    for i in range(n_tmd):
        tmd_phi.append(i*(2*np.pi/n_tmd))
    for (i,t) in enumerate(theta_head1):
        n_phi=int(2*np.pi*r1*np.sin(t)/dpt)+1
        phi=np.linspace(0,2*np.pi,n_phi)
        phi=((phi+np.roll(phi,-1))/2)[0:-1]
        dp_lip = phi[1]-phi[0]
        for (j,p) in enumerate(phi):
                
                    
            x=r1*math.sin(t)*math.cos(p)
            y=r1*math.sin(t)*math.sin(p)
            z=r1*math.cos(t)
            xh.append(x)
            yh.append(y)
            zh.append(z)
    
        for i in range(n_tmd*ec_lef):
            phi = tmd_phi[i]
            x = r_tmd*np.cos(phi)
            y = r_tmd*np.sin(phi)
            z = z_tmd
            xtmd.append(x)
            ytmd.append(y)
            ztmd.append(z)
    
    xh=np.array(xh)
    yh=np.array(yh)
    zh=np.array(zh)
    xtmd=np.array(xtmd)
    ytmd=np.array(ytmd)
    ztmd=np.array(ztmd)
    return(xh,yh,zh,xtmd,ytmd,ztmd)

##### Define a planer bilayer in system
# NOT USED
def planar_mem(box_dim,z_mid,rho,beads_per_lipid,r_tail,r_head):
    hmem = 2 * ((beads_per_lipid - 1) * 2 * r_tail + 2 * r_head) - 2*r_head
    dxy=np.sqrt(1/rho) # distance along between lipids in x and y direction at given density
    bx=box_dim[0]
    by=box_dim[1]
    bz=box_dim[2]
    alipid=1/rho
    n_x=int(bx/dxy)
    n_y=int(by/dxy)
    pos_x = np.zeros(beads_per_lipid *40* n_x * n_y)
    pos_y = np.zeros(beads_per_lipid *40* n_x * n_y)
    pos_z = np.zeros(beads_per_lipid *40* n_x * n_y)
    bx=n_x*dxy
    by=n_y*dxy
    x_head=np.array(2*[(i+0.5)*dxy for i in range(n_x) for j in range(n_y)])-bx/2.0
    y_head=np.array(2*[(j+0.5)*dxy for i in range(n_x) for j in range(n_y)])-by/2.0
    zh = [hmem/2.0+z_mid for k in range(int(len(x_head)/2))]
    z_head=zh+[-hmem/2.0+z_mid for k in range(int(len(x_head)/2))]
    N_lipid=len(x_head)
    N_lipid_spec = int(bx*by/alipid)*2
    if N_lipid < N_lipid_spec:
        raise Exception('Number of lipid mismatch %d placed, %d needed' % (N_lipid,N_lipid_spec))
        print('DEBUG B1: Number of lipid check if matches: %d regularly placed, %d needed' % (N_lipid,N_lipid_spec))
    n1 = [0 for k in x_head]
    n2 = [0 for k in x_head]
    n3 = [-1 for k in range(int(len(x_head)/2))]
    n3 = n3 + [1 for k in range(int(len(x_head)/2))]
    zh = np.array(zh)
    n1 = np.array(n1)
    n2 = np.array(n2)
    n3 = np.array(n3)
    indx=0
    for k in range(len(x_head)):
        x,y,z = rotate_lipid(x_head[k],y_head[k],z_head[k],np.array([n1[k],n2[k],n3[k]]),r_head,r_tail,beads_per_lipid)
        for m in np.arange(len(x)):
            pos_x[indx] = x[m]
            pos_y[indx] = y[m]
            pos_z[indx] = z[m]
            indx = indx + 1
    pos_x = pos_x[0:indx]
    pos_y = pos_y[0:indx]
    pos_z = pos_z[0:indx]
    print('DEBUG B2: Number of beads check if matches: %d needed times %d beads per lipid = %d beads' % (N_lipid_spec,beads_per_lipid,indx))
    return(pos_x,pos_y,pos_z,bx,by,bz)

def rotate_lipid(xh,yh,zh,nrm,r_head,r_tail,beads_per_lipid):

    '''
    Given head coordinates xh,yh,zh
    Orientation of lipid given by unit vector nrm = [nrm_x,nrm_y,nrm_z]
    head radius and tail radius and beads per lipid given by resp. variables
    return x,y,z of the lipid beads
    '''

    # set up the arrays and initialize first beads to head position
    xl = np.zeros(beads_per_lipid)
    yl = np.zeros(beads_per_lipid)
    zl = np.zeros(beads_per_lipid)

    xl[0] = xh # store head beads x coor as the first element in xl array
    yl[0] = yh # store head beads y coor as the first element in yl array
    zl[0] = zh # store head beads z coor as the first element in zl array
    # calculate positions of the other beads-->(tail beads)
    for i in np.arange(beads_per_lipid - 1): # equaverlent to range(beads_per_lipid - 1)
        i = i + 1
        if i == 1:
            xl[i] = xl[i-1] + (r_head + r_tail) * nrm[0] # x coor for first tail beads, store as second element in xl
            yl[i] = yl[i-1] + (r_head + r_tail) * nrm[1] # y coor for first tail beads, store as second element in yl
            zl[i] = zl[i-1] + (r_head + r_tail) * nrm[2] # z coor for first tail beads, store as second element in zl
        else:
            xl[i] = xl[i-1] + (2 * r_tail) * nrm[0]
            yl[i] = yl[i-1] + (2 * r_tail) * nrm[1]
            zl[i] = zl[i-1] + (2 * r_tail) * nrm[2]
    return xl,yl,zl

######## Define pair potentials ########
# ATTRACTIVE COSINE POTENTIAL
def cos_potential(r, rmin, rmax, epsilon, rc, wc):
    V = 0
    F = 0
    if r < rc :
        V = -epsilon
        F = 0
    elif r > rc and r < rc + wc :
        V = -epsilon * (math.cos(math.pi * (r - rc) / (2*wc)) ** 2)
        F = -epsilon * math.pi / wc * math.cos(math.pi * (r - rc) / (2*wc)) * math.sin(math.pi * (r - rc) / (2*wc))

    return (V, F)

######## Place fusogen center beads ########
def place_rod_center(linker_tmd_pair,rod_bd_dia,n_bds_per_rod,rod_len):
    n_rods = len(linker_tmd_pair)
    rod_orien_i = np.array([1,0,0])
    rod_x=[]
    rod_y=[]
    rod_z=[]
    rod_orien=[]
    for i in range(n_rods):
        tmd_x1 = linker_tmd_pair[i][0][0]
        tmd_x2 = linker_tmd_pair[i][1][0]
        tmd_y1 = linker_tmd_pair[i][0][1]
        tmd_y2 = linker_tmd_pair[i][1][1]
        tmd_z1 = linker_tmd_pair[i][0][2]
        tmd_z2 = linker_tmd_pair[i][1][2]
        x =  np.average([tmd_x1,tmd_x2])
        y = np.average([tmd_y1,tmd_y2])
        z = np.average([tmd_z1,tmd_z2])
        r = np.sqrt((x**2)+(y**2))
        r_rod_c = r+rod_len*0.5#((n_bds_per_rod)/2)*rod_bd_dia
        x_c = r_rod_c*(x/r)
        y_c = r_rod_c*(y/r)
        nrm=np.array([x_c,y_c,0])
        nrm=nrm/(np.sqrt((nrm[0]**2)+(nrm[1]**2)+(nrm[2]**2)))
        rot_axis = np.cross(nrm,rod_orien_i)
        if sum(rot_axis) == 0:
            #orientation = [1.0,0.0,0.0,0.0]
            nrm_dot_rod_orien = np.dot(nrm,rod_orien_i)
            nrm_dot_rod_orien = np.dot(nrm,rod_orien_i)
            rot_ang = -np.arccos(nrm_dot_rod_orien)
            rod_x.append(x_c)
            rod_y.append(y_c)
            rod_z.append(z)
            orientation=[np.cos(rot_ang/2),rot_axis[0]*np.sin(rot_ang/2),rot_axis[1]*np.sin(rot_ang/2),rot_axis[2]*np.sin(rot_ang/2)]
        else:
            rot_axis=rot_axis/(np.sqrt((rot_axis[0]**2)+(rot_axis[1]**2)+(rot_axis[2]**2)))
            nrm_dot_rod_orien = np.dot(nrm,rod_orien_i)
            nrm_dot_rod_orien = np.dot(nrm,rod_orien_i)
            rot_ang = -np.arccos(nrm_dot_rod_orien)
            rod_x.append(x_c)
            rod_y.append(y_c)
            rod_z.append(z)
            orientation=[np.cos(rot_ang/2),rot_axis[0]*np.sin(rot_ang/2),rot_axis[1]*np.sin(rot_ang/2),rot_axis[2]*np.sin(rot_ang/2)]
        rod_orien.append(orientation)
    return(np.array(rod_x),np.array(rod_y),np.array(rod_z),np.array(rod_orien))

def linker_pot(r,rmin,rmax,k,r1):
    V = 0
    F = 0
    if r < r1 :
        V = 0.5*k*(r**2)
        F = -k*r
    elif r >= r1 :
        V = k*r1*r+(0.5*k*(r1**2)-k*r1*r1)
        F = -k*r1
    return (V,F)

def linker_pot_real_snare(r,rmin,rmax,N_unzip,lp,linker):
    lp = lp/0.88 
    a = 0.365/0.88 # N_unzip is the number of uncomplexed residues
    rr = r #- (1+0.66)/0.88/2
    if linker == 0:
        V = 0
        F = 0
    else:
        if rr > 0:
            V = a*N_unzip/(4*lp*(1-rr/(a*N_unzip)))*(3*(rr/(a*N_unzip))**2-2*(rr/(a*N_unzip))**3)/0.6
            F = -1/lp*(1/(4*(1-rr/(a*N_unzip))**2)+rr/(a*N_unzip)-1/4)/0.6
        else:
            V = 0
            F = 0
    return (V,F)

# NOT USED
def elect_pot(r,rmin,rmax,r_bead1,r_bead2,charge1,charge2,wcharge):
    lambda_d = 0.8/0.88 # unit in sigma
    if ((r_bead1 == 0) or (r_bead2 == 0)):
        V = 0
        F = 0
    else:
        r_eff = r_bead1*r_bead2/(r_bead1+r_bead2) # unit in sigma
        V = 0
        F = 0
        if wcharge!=0:
            e = 1.6*10**(-19)
            Psi1 = 53.4*math.asinh(charge1*e/(4*math.pi*(r_bead1*0.88*10**(-9))**2)/(0.116*0.15**0.5)) # unit of charge density is in C/m^2
            Psi2 = 53.4*math.asinh(charge2*e/(4*math.pi*(r_bead2*0.88*10**(-9))**2)/(0.116*0.15**0.5))
            Z1 = 9.38*10**(1)*math.tanh(Psi1/107)**2*0.88/4.1/0.6 # unit in epsilon/sigma
            Z2 = 9.38*10**(1)*math.tanh(Psi2/107)**2*0.88/4.1/0.6 # unit in epsilon/sigma
            if (r-r_bead1-r_bead2)/lambda_d>0:
                if r-r_bead1-r_bead2 >=0:
                    if charge1*charge2 < 0:
                        V = -r_eff*math.sqrt(abs(Z1*Z2))*math.exp(-(r-r_bead1-r_bead2)/lambda_d)
                        F = -r_eff/lambda_d*math.sqrt(abs(Z1*Z2))*math.exp(-(r-r_bead1-r_bead2)/lambda_d)
                    else:
                        V = r_eff*math.sqrt(abs(Z1*Z2))*math.exp(-(r-r_bead1-r_bead2)/lambda_d)
                        F = r_eff/lambda_d*math.sqrt(abs(Z1*Z2))*math.exp(-(r-r_bead1-r_bead2)/lambda_d)
                else:
                    if charge1*charge2 < 0:
                        V = -r_eff*math.sqrt(abs(Z1*Z2))
                        F = -r_eff/lambda_d*math.sqrt(abs(Z1*Z2))*0
                    else:
                        V = r_eff*math.sqrt(abs(Z1*Z2))
                        F = r_eff/lambda_d*math.sqrt(abs(Z1*Z2))*0
    return (V,F) 

# NOT USED
def elect_pot_screened(r,rmin,rmax,r_bead1,r_bead2,charge1,charge2,wcharge):
    lambda_d = 0.8/0.88 # unit in sigma
    if ((r_bead1 == 0) or (r_bead2 == 0)):
        V = 0
        F = 0
    else:
        r_eff = r_bead1*r_bead2/(r_bead1+r_bead2) # unit in sigma
        V = 0
        F = 0
        if wcharge!=0:
            e = 1.6*10**(-19)
            Psi1 = 53.4*math.asinh(charge1*e/(4*math.pi*(r_bead1*0.88*10**(-9))**2)/(0.116*0.15**0.5)) # unit of charge density is in C/m^2
            Psi2 = 53.4*math.asinh(charge2*e/(4*math.pi*(r_bead2*0.88*10**(-9))**2)/(0.116*0.15**0.5))
            Z1 = 9.38*10**(1)*math.tanh(Psi1/107)**2*0.88/4.1/0.6 # unit in epsilon/sigma
            Z2 = 9.38*10**(1)*math.tanh(Psi2/107)**2*0.88/4.1/0.6 # unit in epsilon/sigma
            if (r-r_bead1-r_bead2)/lambda_d>0:
                if r-r_bead1-r_bead2 >=0:
                    if charge1*charge2 < 0:
                        V = -r_eff*math.sqrt(abs(Z1*Z2))*math.exp(-(r-r_bead1-r_bead2)/lambda_d)
                        F = -r_eff/lambda_d*math.sqrt(abs(Z1*Z2))*math.exp(-(r-r_bead1-r_bead2)/lambda_d)
                    else:
                        V = r_eff*math.sqrt(abs(Z1*Z2))*math.exp(-(r-r_bead1-r_bead2)/lambda_d)
                        F = r_eff/lambda_d*math.sqrt(abs(Z1*Z2))*math.exp(-(r-r_bead1-r_bead2)/lambda_d)
                else:
                    if charge1*charge2 < 0:
                        V = -r_eff*math.sqrt(abs(Z1*Z2))
                        F = -r_eff/lambda_d*math.sqrt(abs(Z1*Z2))*0
                    else:
                        V = r_eff*math.sqrt(abs(Z1*Z2))
                        F = r_eff/lambda_d*math.sqrt(abs(Z1*Z2))*0
    return (V,F) 

############################################################
################## SET IC SECTION ##########################
############################################################

##### scan parameters #####
### num of TMD hydrophobic beads

n_hypho_bds_tmd = 2

### strength of tmd-lip attr

const = 3

### N staple strength

const_n=2

### C staple strength

const_c=2

### constant force magnitude

f_star_grp = [18,20,22,25,30,35]

### TMD bead size

tmd_bd_size_nm = 1
tmd_bd_r_nm = 0.5*tmd_bd_size_nm #tmd bead size in nm
tmd_bd_r = tmd_bd_r_nm/0.88

### Number of SNARE

n_rod_grp = [12,9,7,6,5,4,3]
n_rod = n_rod_grp[int(sys.argv[5])]

### Length of SNARE

rod_len_nm=10

### Membrane Tension

tension = 0.05 #pN/nm

### run number

run_i = int(sys.argv[1])+0
run = run_i

### number of unzippered residues and beads

N_unzip_grp = [7,10,14,18]
N_unzip = N_unzip_grp[int(sys.argv[6])]
N_bd_unzip = 10//3##N_unzip//3

### Lp setting 

lp_grp = [0.3,0.35,0.4,0.45,0.5,0.55,0.6]
lp = lp_grp[int(sys.argv[7])]

### mean force critrica

mean_force_grp = [33]
mean_force = mean_force_grp[int(sys.argv[1])//1%1]

# tmd beads types and positions, get from example_IC_7HA2_10nm_LR_processed_beads.csv
tmd_hypho_bd_dia = tmd_bd_r*2
tmd_hyphi_bd_dia = tmd_bd_r*2
n_hyphi_bds_tmd = 2
n_bds_tmd = int(n_hypho_bds_tmd + n_hyphi_bds_tmd)
tmd_hypho_len = (n_hypho_bds_tmd-1)*tmd_hypho_bd_dia
tmd_hypho_x = np.array([tmd_hypho_bd_dia,0.5*tmd_hypho_bd_dia,-0.5*tmd_hypho_bd_dia])
#tmd_hypho_x = np.linspace(0.5*tmd_hypho_len,-0.5*tmd_hypho_len,n_hypho_bds_tmd)
tmd_hyphi_x = np.array([tmd_hypho_x[1]+0.5*tmd_hypho_bd_dia+0.5*tmd_hyphi_bd_dia,tmd_hypho_x[-1]-0.5*tmd_hypho_bd_dia-0.5*tmd_hyphi_bd_dia])
tmd_x = np.append(tmd_hyphi_x,tmd_hypho_x)
tmd_y = np.zeros(n_bds_tmd+1)
tmd_z = np.zeros(n_bds_tmd+1)
tmd_types = ['N' for i in range(n_hyphi_bds_tmd-1)]+['C']+['LT']+['E' for j in range(n_hypho_bds_tmd)]
print(tmd_types)

### Fusogen beads positions and types, based on the fusogen setup note
### Place the rod beads along x-axis so that x-axis would be the body axis and centered at COM
### sys.argv[2] is 0 for SNARE runs
if int(sys.argv[2])==0:
    n_bds_per_bundle = [14+3-N_bd_unzip,16,20,19+3-N_bd_unzip,14+3-N_bd_unzip,16,20,19+3-N_bd_unzip] # add N-terminal beads for VAMP and Syntaxin 

# The first bead of the bindles corresponding to bead 20 represents layer +8 of the SNAREpins
    rod_type_vamp = ['0J']*(3-N_bd_unzip)+['1J','n1J','0J','0J','n1J','n2J','0J','0J','0J','1J','n2J','n1J','0J','2J']
    rod_type_sn2 = ['n1J','1J','n2J','2J','0J','0J','0J','n1J','n1J','n1J','1J','0J','0J','n2J','n2J','0J']
    rod_type_sn1 = ['n1J','0J','0J','n1J','0J','n1J','n2J','n2J','0J','1J','0J','n2J','1J','0J','n2J','n1J','2J','n2J','1J','0J']
    rod_type_stx = ['0J']*(3-N_bd_unzip)+['2J','n1J','0J','n1J','n1J','0J','n2J','n1J','0J','n1J','n1J','0J','0J','n1J','0J','n1J','2J','0J','0J']
    rod_type_vamp_2 = ['LL']*(3-N_bd_unzip)+['VA' for i in range(14)]
    rod_type_sn2_2 = ['SN2' for i in range(16)]
    rod_type_sn1_2 = ['SN1' for i in range(20)]
    rod_type_stx_2 = ['LL']*(3-N_bd_unzip)+['STX' for i in range(19)]
    c_term_types =['LJ','LJ']
    c_term_types_2 =['L','L']
    c_term_types_i =['LJI','LJI']
    c_term_types_2_i =['LI','LI']
    rod_x = []
    rod_y = []
    rod_z = []
    rod_bd_dia_nm = 0.66
    rod_bd_dia = rod_bd_dia_nm/0.88
    rod_len = rod_len_nm/0.88
    rod_length = rod_len - rod_bd_dia

    for ii in range(8):
        n_bds_per_cyl = n_bds_per_bundle[ii]
        if ii in [0,3,4,7]:
            cyl_x = [2*rod_bd_dia/2*i for i in range(n_bds_per_cyl)]
            cyl_y = [np.sqrt(2)*1.33/0.88/2*math.cos((-14*i-225-90*ii)/180*math.pi) for i in range(n_bds_per_cyl)]
            cyl_z = [np.sqrt(2)*1.33/0.88/2*math.sin((-14*i-225-90*ii)/180*math.pi) for i in range(n_bds_per_cyl)]
        else:
            cyl_x = [2*rod_bd_dia/2*i for i in range(3-N_bd_unzip,n_bds_per_cyl+3-N_bd_unzip)]
            cyl_y = [np.sqrt(2)*1.33/0.88/2*math.cos((-14*i-225-90*ii)/180*math.pi) for i in range(3-N_bd_unzip,n_bds_per_cyl+3-N_bd_unzip)]
            cyl_z = [np.sqrt(2)*1.33/0.88/2*math.sin((-14*i-225-90*ii)/180*math.pi) for i in range(3-N_bd_unzip,n_bds_per_cyl+3-N_bd_unzip)]
        rod_x = rod_x + cyl_x
        rod_y = rod_y + cyl_y
        rod_z = rod_z + cyl_z
    
    # vamp is ii = 0, syntaxin is ii = 3
    ii = 0
    i = -1
    c_term_x_vamp = 2*rod_bd_dia/2*i
    c_term_y_vamp = np.sqrt(2)*1.33/0.88/2*math.cos((-14*i-225-90*ii)/180*math.pi)
    c_term_z_vamp = np.sqrt(2)*1.33/0.88/2*math.sin((-14*i-225-90*ii)/180*math.pi)
    
    ii = 3
    i = -1
    c_term_x_stx = 2*rod_bd_dia/2*i
    c_term_y_stx = np.sqrt(2)*1.33/0.88/2*math.cos((-14*i-225-90*ii)/180*math.pi)
    c_term_z_stx = np.sqrt(2)*1.33/0.88/2*math.sin((-14*i-225-90*ii)/180*math.pi)

    c_term_x = [c_term_x_vamp,c_term_x_stx]
    c_term_y = [c_term_y_vamp,c_term_y_stx]
    c_term_z = [c_term_z_vamp,c_term_z_stx]
    
    rod_x = rod_x + c_term_x + c_term_x
    rod_y = rod_y + c_term_y + c_term_y
    rod_z = rod_z + c_term_z + c_term_z
    rod_x = np.array(rod_x)
    rod_y = np.array(rod_y)
    rod_z = np.array(rod_z)
    rod_x = rod_x-(np.max(rod_x)-np.min(rod_x))/2
    rod_types=rod_type_vamp+rod_type_sn2+rod_type_sn1+rod_type_stx+rod_type_vamp_2+rod_type_sn2_2+rod_type_sn1_2+rod_type_stx_2+c_term_types+c_term_types_2
    rod_types_i=rod_type_vamp+rod_type_sn2+rod_type_sn1+rod_type_stx+rod_type_vamp_2+rod_type_sn2_2+rod_type_sn1_2+rod_type_stx_2+c_term_types_i+c_term_types_2_i
    n_bds_per_rod = len(rod_types)
else:
    rod_bd_dia_nm = 2#2.5
    rod_bd_dia = rod_bd_dia_nm/0.88
    rod_bd_dia_fix = 2/0.88
    c_term_dia = (np.sqrt(2)-1)*rod_bd_dia_fix/(np.sqrt(2)+1)

    n_bds_per_cyl = rod_len_nm-1#9

    rod_len = rod_len_nm/0.88
    rod_length = rod_len - rod_bd_dia #(n_bds_per_cyl-1-4) * rod_bd_dia
    cyl_x=np.linspace(-0.5*rod_length,0.5*rod_length,n_bds_per_cyl) #np.linspace(-(10/0.88)*0.5,(10/0.88)*0.5,n_bds_per_cyl)#np.linspace(-0.5*rod_length,0.5*rod_length,n_bds_per_cyl)
    cyl_y=np.zeros(n_bds_per_cyl)
    cyl_z=np.zeros(n_bds_per_cyl)

    c_term_x_i=cyl_x[0]-np.sqrt(((0.5*rod_bd_dia+tmd_bd_r)**2)-(tmd_bd_r**2))
    c_term_x=np.array([c_term_x_i,c_term_x_i])
    c_term_y = np.array([0,0])

    c_term_z = np.array([tmd_bd_r,-tmd_bd_r])
    rod_x = np.concatenate((c_term_x,cyl_x))
    rod_y = np.concatenate((c_term_y,cyl_y))
    rod_z = np.concatenate((c_term_z,cyl_z))
    cyl_types=['R' for i in range(n_bds_per_cyl)]
    c_term_types =['L','L']
    rod_types=c_term_types+cyl_types
    n_bds_per_rod = len(rod_types)

# append each rod bead into an array
rod_pos=[]
for i in range(len(rod_x)):
    x=rod_x[i]
    y=rod_y[i]
    z=rod_z[i]
    rod_pos.append([x,y,z])
print('rod beads position:')
print(rod_pos)

# append each tmd bead into an array
tmd_pos=[]
for i in range(len(tmd_x)):
    x=tmd_x[i]
    y=tmd_y[i]
    z=tmd_z[i]
    tmd_pos.append([x,y,z])

print('tmd beads position:')
print(tmd_pos)
NTMDBeads=len(tmd_types)

# Set up moment of inertia of TMDs and SNAREpins
m_tmd = 1 
l_tmd_x = n_bds_tmd*tmd_bd_r*2
l_tmd_y = tmd_bd_r
l_tmd_z = tmd_bd_r
i_tmd_x = 0.5*m_tmd*((0.5*(l_tmd_y+l_tmd_z))**2)
i_tmd_y = (1/12)*m_tmd*(3*((0.5*(l_tmd_y+l_tmd_z))**2)+l_tmd_x**2)
i_tmd_z = i_tmd_y
tmd_len_nm = l_tmd_x*0.88

m_rod = 1
if int(sys.argv[2])==0:
    l_rod_x = max(rod_x)-min(rod_x)+0.66/0.88#rod_bd_dia #max(rod_x)-min(rod_x)
else:
    l_rod_x = max(rod_x)-min(rod_x)+2/0.88
l_rod_y = 0.5*2/0.88#rod_bd_dia#max(rod_y)-min(rod_y)
l_rod_z = 0.5*2/0.88#rod_bd_dia#max(rod_z)-min(rod_z)
i_rod_x = 0.5*m_rod*((0.5*(l_rod_y+l_rod_z))**2)
i_rod_y = (1/12)*m_rod*(3*((0.5*(l_rod_y+l_rod_z))**2)+l_rod_x**2)
i_rod_z = i_rod_y

# Set up membrane tension
density = (tension-46.9)/(-63.7) # Density of lipids corresponds to specified tension
sigma_bead=1
r_ves=21
r_head = 0.95 * sigma_bead / 2.0 # Head Radius
r_tail = 1 * sigma_bead / 2.0 # Tail Radius
beads_per_lipid = 4 # Number of beads per lipid

r_eci_arr = [8.5,8,7.5,7,6.5,6,5.5,5,4.5,4,3.5,3,2.5,2]
r_eci = r_eci_arr[(int(sys.argv[1])-1)]
if r_eci >= 7.0:
    delta_h = 1.0
elif r_eci >= 5:
    delta_h = 1.5
elif r_eci >= 4:
    delta_h = 1.8
elif r_eci >= 3.2:
    delta_h = 1.85
elif r_eci >= 2:
    delta_h = 2.1

hmem = 2 * ((beads_per_lipid - 1) * 2 * r_tail + 2 * r_head) - 2*r_head
z_c_1=r_ves+hmem+delta_h
z_c_2=-(r_ves+hmem+delta_h)

n_tmd=n_rod
tmd_theta=np.pi-np.arcsin(r_eci/(r_ves+hmem))
n_gas_i=int(tension*8*np.pi*(r_ves**2)/3/4)+1
n_gas = n_gas_i
r_gas=1*sigma_bead
mass_lip=1
mass_gas=mass_lip

###
### INITIAL PARAMETERS
###
sigma_bead = float(1.0) # TAIL BEAD SIZE (DIAMETER)
alipid = float(1/density) # AREA PER LIPID
rhead = 0.95 * sigma_bead / 2.0 # Head Radius
rtail = 1 * sigma_bead / 2.0 # Tail Radius
bx=100 # box size in x
by=100 # box size in y
bz=140 # box size in z
# Get lipid beads locations and tmds locations
pos_x,pos_y,pos_z,tmd_c_x,tmd_c_y,tmd_c_z,tmd_orien=vesicle(r_ves,z_c_1,density,beads_per_lipid,r_tail,r_head,n_tmd,tmd_theta) # top vesicel
pos_x_2,pos_y_2,pos_z_2,tmd_c_x_2,tmd_c_y_2,tmd_c_z_2,tmd_orien_2=vesicle(r_ves,z_c_2,density,beads_per_lipid,r_tail,r_head,n_tmd,np.pi-tmd_theta) # bottom vesicel
pos_x = np.concatenate((pos_x,pos_x_2))
pos_y = np.concatenate((pos_y,pos_y_2))
pos_z = np.concatenate((pos_z,pos_z_2))
tmd_c_x = np.concatenate((tmd_c_x,tmd_c_x_2))
tmd_c_y = np.concatenate((tmd_c_y,tmd_c_y_2))
tmd_c_z = np.concatenate((tmd_c_z,tmd_c_z_2))
tmd_orien.extend(tmd_orien_2)
# create a snapshot for the central ghost beads
b_type_arr=['bond1','bond2','bond3','bond4','linker','linker2']
if int(sys.argv[2])!=0:
    p_type_arr=['H','T','G','D','N','C','E','A','L','R']
else:
    p_type_arr=['H','T','G','D','N','N2','C','E','A','n2J','n1J','0J','1J','2J','VA','STX','SN1','SN2','LL','3J','L','LJ','LT','D2','A2','LI','LJI'] # D A are the center beads
print(len(tmd_types))
snap = hoomd.data.make_snapshot(N=len(tmd_c_x)+n_rod, box=hoomd.data.boxdim(Lx=bx,Ly=by,Lz=bz,dimensions=3,volume=None), particle_types=p_type_arr,bond_types=b_type_arr)

# place beads
ang = 0
threeFlag = 0

tmd_pos_r = []
# Place TMDs center beads
for k in range(len(tmd_c_x)):
    snap.particles.position[k] = [tmd_c_x[k],tmd_c_y[k],tmd_c_z[k]]
    snap.particles.orientation[k] = tmd_orien[k]
    if k == len(tmd_c_x)-1 or k == int(len(tmd_c_x)/2-1):
        tmd_c_index=p_type_arr.index('D2')
    else:
        tmd_c_index=p_type_arr.index('D')
    snap.particles.typeid[k] = tmd_c_index
    snap.particles.moment_inertia[k] = [i_tmd_x,i_tmd_y,i_tmd_z]
    snap.particles.mass[k] = 1  #1*NTMDBeads # need to change
    r_tmd = np.sqrt((tmd_c_x[k]**2)+(tmd_c_y[k]**2))
    tmd_pos_r.append(r_tmd)
tmd_pos_r = np.array(tmd_pos_r)
tmd_pos_r_avg = np.mean(tmd_pos_r)
print('mean tmd center radius:')
print(tmd_pos_r_avg)

# Place rods center beads
for k in range(n_rod):
    index = k+len(tmd_c_x)
    if k < n_rod-1:
        rod_c_index=p_type_arr.index('A')
    else:
        rod_c_index=p_type_arr.index('A2')
    snap.particles.typeid[index] = rod_c_index
    snap.particles.moment_inertia[index] = [i_rod_x,i_rod_y,i_rod_z]
    snap.particles.mass[index] = 1   #1*n_bds_per_cyl

system = hoomd.init.read_snapshot(snap)
print(snap.particles.N)

# Add rigid body constituent particles

# Place TMD constitutuent beads
rigid = hoomd.md.constrain.rigid()
rigid.set_param('D', types=tmd_types, positions=tmd_pos)
rigid.set_param('D2', types=tmd_types, positions=tmd_pos)
rigid.create_bodies()
n_tmd_beads = len(system.particles)
print("Number of particles after create %.0f TMDs is %.1f" %(n_tmd,n_tmd_beads))
n_tmd_c = n_tmd*2
n_rod_c = n_rod
linker_tmds_pos=[]
snap=system.take_snapshot()
for i in range(n_tmd):
    tmd_tag_1 = n_tmd_c+n_rod_c+i*NTMDBeads
    tmd_tag_2 = tmd_tag_1+n_tmd*NTMDBeads
    tmd_1_pos=snap.particles.position[tmd_tag_1]
    tmd_2_pos=snap.particles.position[tmd_tag_2]
    linker_tmds_pos.append([tmd_1_pos,tmd_2_pos])
print('linker_tmd pos:')
print(linker_tmds_pos)
print(l_rod_x)
rod_c_x,rod_c_y,rod_c_z,rod_orien=place_rod_center(linker_tmds_pos,rod_bd_dia,n_bds_per_cyl,l_rod_x)
index = 0
for p in system.particles:
    if p.type == 'A' or p.type == 'A2':
        p.position=[rod_c_x[index],rod_c_y[index],rod_c_z[index]]
        print([rod_c_x[index],rod_c_y[index],rod_c_z[index]])
        p.orientation = rod_orien[index]
        p.moment_inertia=[i_rod_x,i_rod_y,i_rod_z]
        print('mass of rod is:')
        print(p.mass)
        index+=1
    elif p.type =='D' or p.type == 'D2':
        print('mass of tmd is:')
        print(p.mass)
rigid.set_param('A', types=rod_types, positions=rod_pos)
rigid.set_param('A2', types=rod_types_i, positions=rod_pos)
rigid.create_bodies()
rigid.validate_bodies()
n_rigid_beads = len(system.particles)
#f_name = 'scan_rod_len_%.1f_tmd_n_attr_%.1f_tmd_c_attr_%.1f_tmd_bd_size_%.2f_n_phobic_bds_%.0f_n_philic_bds_%.0f_nrod_%.0f_rod_r_%.2f_eliptmd_%.2f_f_star_%.1f_dves_%.1f_tension_%.2f_run_%.0f' % (rod_len_nm,const_n,const_c,tmd_bd_r,n_hypho_bds_tmd,n_hyphi_bds_tmd,n_rod,rod_bd_dia/2,const,f_star_pN,2*(r_ves+hmem),tension,run)
#hoomd.dump.gsd(f_name+'_ini.gsd',period=None,group=hoomd.group.all(),overwrite=True)
#exit()
snap=system.take_snapshot()
print(snap.particles.N)
# Add linkers
linker_length=[]
tmd_d=[]
for i in range(n_rod):
    tmd_tag_1 = n_tmd_c+n_rod_c+i*NTMDBeads+2
    tmd_tag_2 = tmd_tag_1 + n_tmd*NTMDBeads+0
    if int(sys.argv[2])==0:
        rod_tag_1 = n_tmd_beads + (i+1)*n_bds_per_rod-4 # The first VAMP bead
        rod_tag_2 = rod_tag_1 + 1 # The first syntaxin bead
    else:
        rod_tag_1 = n_tmd_beads + i*n_bds_per_rod
        rod_tag_2 = rod_tag_1 + 1
    if i < n_rod - 1:
        system.bonds.add('linker',rod_tag_1,tmd_tag_1)
        system.bonds.add('linker',rod_tag_2,tmd_tag_2)
    else:
        system.bonds.add('linker2',rod_tag_1,tmd_tag_1)
        system.bonds.add('linker2',rod_tag_2,tmd_tag_2)
    
    tmd_1_pos = np.array(system.particles[tmd_tag_1].position)
    tmd_2_pos = np.array(system.particles[tmd_tag_2].position)
    rod_pos1 = np.array(system.particles[rod_tag_1].position)
    rod_pos2 = np.array(system.particles[rod_tag_2].position)
    
    length_1 = np.sqrt(((rod_pos1[0]-tmd_1_pos[0])**2)+((rod_pos1[1]-tmd_1_pos[1])**2)+((rod_pos1[2]-tmd_1_pos[2])**2))
    length_2 = np.sqrt(((rod_pos2[0]-tmd_2_pos[0])**2)+((rod_pos2[1]-tmd_2_pos[1])**2)+((rod_pos2[2]-tmd_2_pos[2])**2))

    d = np.sqrt(((tmd_1_pos[0]-tmd_2_pos[0])**2)+((tmd_1_pos[1]-tmd_2_pos[1])**2)+((tmd_1_pos[2]-tmd_2_pos[2])**2))
    linker_length.append([length_1,length_2])
    tmd_d.append(d)
print('lengths of linkers are:')
print(linker_length)
#f_name = 'scan_rod_len_%.1f_tmd_n_attr_%.1f_tmd_c_attr_%.1f_tmd_bd_size_%.2f_n_phobic_bds_%.0f_n_philic_bds_%.0f_nrod_%.0f_rod_r_%.2f_eliptmd_%.2f_f_star_%.1f_dves_%.1f_tension_%.2f_run_%.0f' % (rod_len_nm,const_n,const_c,tmd_bd_r,n_hypho_bds_tmd,n_hyphi_bds_tmd,n_rod,rod_bd_dia/2,const,f_star_pN,2*(r_ves+hmem),tension,run)

#print(f_name)
#hoomd.dump.gsd(f_name+'_ini.gsd',period=None,group=hoomd.group.all(),overwrite=True)
#exit()
### Check lipid-tmd beads overlapping and remove overlapping lipids
tmd_rod_snap=system.take_snapshot()
tmd_rod_pos=tmd_rod_snap.particles.position
tmd_bd_types=['N','C','E','LT']

tmd_pos_arr=[]

for (i,pos) in enumerate(tmd_rod_pos):
    bd_typeid=tmd_rod_snap.particles.typeid[i]
    bd_type=p_type_arr[bd_typeid]
    if bd_type in tmd_bd_types:
        tmd_pos_arr.append(pos)
    

# search for tmd-overlapping lipids
lip_i_remov_arr=[]
for (i,tmd_pos) in enumerate(tmd_pos_arr):
    tmd_x=tmd_pos[0]
    tmd_y=tmd_pos[1]
    tmd_z=tmd_pos[2]
    for (j,lip_x) in enumerate(pos_x):
        lip_y=pos_y[j]
        lip_z=pos_z[j]
        lip_tmd_r_sq=((tmd_x-lip_x)**2)+((tmd_y-lip_y)**2)+((tmd_z-lip_z)**2)
        if lip_tmd_r_sq<=(0.5*sigma_bead+tmd_bd_r)**2:
            lip_i=j//beads_per_lipid
            lip_i_remov_arr.append(lip_i)

# remove tmd-overlapping lipid beads positions
lipid_pos_x=[]
lipid_pos_y=[]
lipid_pos_z=[]
for (i,lip_x) in enumerate(pos_x):
    if i//beads_per_lipid in lip_i_remov_arr:
        pass
    else:
        lip_y=pos_y[i]
        lip_z=pos_z[i]
        lipid_pos_x.append(lip_x)
        lipid_pos_y.append(lip_y)
        lipid_pos_z.append(lip_z)
'''
for p in system.particles:
    if p.type == 'N':
        p_pos = p.position
        p_tag = p.tag
        if p_pos[2]<0:            
            system.particles[p_tag].type='N2'
'''

### Place lipid beads
n_lipid_bds = len(lipid_pos_x)
print('number of lipid beads: %.1f'%n_lipid_bds)
for k in range(n_lipid_bds):
    p_pos = [lipid_pos_x[k],lipid_pos_y[k],lipid_pos_z[k]]
    p_tag=k+n_rigid_beads
    if k%beads_per_lipid==0:
        system.particles.add('H')
        system.particles[p_tag].position=p_pos
    else:
        system.particles.add('T')
        system.particles[p_tag].position=p_pos

# add lipid bonds
n_tmd_bds=n_tmd*NTMDBeads
for p in system.particles:
    if p.type=='H':
        bd_tag=p.tag
        system.bonds.add('bond1',bd_tag,bd_tag+1)
        system.bonds.add('bond2',bd_tag+1,bd_tag+2)
        system.bonds.add('bond3',bd_tag+2,bd_tag+3)
        system.bonds.add('bond4',bd_tag,bd_tag+3)
for i in range(n_gas):
    gas_tag_1=n_rigid_beads+n_lipid_bds+i*2
    system.particles.add('G')
    system.particles[gas_tag_1].position=[0,0,z_c_1]
    system.particles[gas_tag_1].mass=mass_gas
    system.particles.add('G')
    gas_tag_2=gas_tag_1+1
    system.particles[gas_tag_2].position=[0,0,z_c_2]
    system.particles[gas_tag_2].mass=mass_gas


global r1_linker_nm
    
r1_linker_nm = 0.1
r1_linker = r1_linker_nm/0.88
#k_linker_pNnm_arr = [5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9] # in pN/nm
#k_linker = k_linker_pNnm_arr[int(sys.argv[6])]*((0.88**2)/0.6/4) # in eplison/sigma^2 unit
#k_linker = (f_star_pN/r1_linker_nm)*((0.88**2)/0.6/4) # k1_linker in eplison/sigma^2 unit
### output initial condition and check
#f_name = 'nrod_%.0f_rod_r_%.2f_eliptmd_%.2f_f_star_%.1f_dves_%.1f_tension_%.2f_run_%.0f_lp_%.2f_Nbdunzip_%.0f_old_linker_%.0f_jin_fusogen_%.0f_with_electrostatic_%.0f_try_2' % (n_rod,rod_bd_dia/2,const,f_star_pN,2*(r_ves+hmem),tension,run,lp,N_bd_unzip,int(sys.argv[3]),int(sys.argv[2]),int(sys.argv[4]))
f_name = 'ves_ves_avg_nostaple_rtmd_%.2f_deltah_%.2f_Nunzip_%.2f_lp_%.2f_nrod_%.0f_rod_r_%.2f_dves_%.1f_tension_%.2f_run_%.0f' % (tmd_pos_r_avg,delta_h,N_unzip,lp,n_rod,rod_bd_dia/2,2*(r_ves+hmem),tension,run)
print(f_name)
hoomd.dump.gsd(f_name+'_ini.gsd',period=None,group=hoomd.group.all(),overwrite=True)
#exit()

if __name__ == '__main__':

    # Simulation parameters
    tau = 1 # MD Correlation Time
    epsilon = 1 # System Energy Scale
    particle_mass = (tau**2) * epsilon / sigma_bead # Particle Mass
    rmincut = 0.1 * sigma_bead # Minimum separation between beads
    seed_sim = randint(0,1000)

    # Time step of simulation
    dt = 0.005 * tau 
    dt_eq = dt/100
    # Gamma: the langevin drag parameter
    Gamma = tau ** (-1)
    GammaG = Gamma
    # Temperature of the system
    T = 1.7 * epsilon
    N_ratio = 3 # ratio of number of tail beads to number of head beads
    wc = 1.6 * sigma_bead
    
    # Set the parameters for the Weeks-Chandler-Anderson potential
    # since this is very similar to lennard-jones, we call it
    # 'lennard jones shifted' throughout the code
    b_h = (1/2)*0.95 * sigma_bead
    b_t = (1/2)*sigma_bead
    b_g = r_gas
    b_d = 0
    b_d2 = 0
    #b_f = (1/2)*0.6*sigma_bead/0.88
    b_f = tmd_bd_r#tmd_philic_r
    b_n = tmd_bd_r
    b_n2 = tmd_bd_r
    b_c = tmd_bd_r
    b_e = tmd_bd_r#(1/2)*sigma_bead/0.88
    b_a = 0
    b_a2 = 0
    b_l = 0
    b_lj = 0
    b_li = 0
    b_lji = 0
    b_lt = 0
    b_r = rod_bd_dia/2
    b_n2j= 0.33/0.88
    b_n1j= 0.33/0.88
    b_0j= 0.33/0.88
    b_1j= 0.33/0.88
    b_2j= 0.33/0.88
    b_3j= 0.33/0.88
    charge_n2j = -2
    charge_n1j = -1
    charge_0j = 0
    charge_1j = 1
    charge_2j = 2
    charge_3j = 3
    b_va = 0.33/0.88
    b_stx = 0.33/0.88
    b_sn1 = 0.33/0.88
    b_sn2 = 0.33/0.88
    b_ll = 0.33/0.88
    # set rc & wc for both the cosine and the lj potential
    for k in itertools.combinations_with_replacement(p_type_arr,2):
        b_command="b_%s_%s = b_%s + b_%s" % (k[0].lower(),k[1].lower(),k[0].lower(),k[1].lower())
        rc_command="rc_%s_%s = (2 ** (1/6.0)) * b_%s_%s" % (k[0].lower(),k[1].lower(),k[0].lower(),k[1].lower())
        wc_command = "wc_%s_%s = 1.6 * (b_%s + b_%s)" % (k[0].lower(),k[1].lower(),k[0].lower(),k[1].lower())
        exec(b_command)
        exec(rc_command)
        exec(wc_command)
  
    b_h_h = 0.95*sigma_bead
    b_h_t = b_h_h
    rc_h_h = (2**(1/6))*b_h_h
    rc_h_t = (2**(1/6))*b_h_t

    # set bond constant
    k_bond = 30 * epsilon / (sigma_bead ** 2) ### BOND POTENTIAL (FENE)
    k_bend = k_bond / 3.0 ### BENDING STIFFNESS
    r_inf = 1.5 * sigma_bead ### DIVERGENCE LENGTH

    # CREATE NEIGHBOR LIST
    nl = hoomd.md.nlist.cell()

    ###
    ### POTENTIALS
    ###
    # Attractive Cosine Potential	
    #const_arr = np.linspace(1.5,3,4)
    table = hoomd.md.pair.table(width=100000,nlist = nl,name='tail_attr')
    for k in itertools.combinations_with_replacement(p_type_arr,2):
        if k[0]=='T' and k[1]=='T':
            table.pair_coeff.set('T', 'T', func=cos_potential, rmin=0,rmax=rc_t_t+wc,coeff=dict(epsilon=epsilon, rc=rc_t_t, wc=wc))
        elif k[0]=='T' and k[1]=='E':
            table.pair_coeff.set('T', 'E', func=cos_potential, rmin=0,rmax=rc_t_e+wc_t_e,coeff=dict(epsilon=const*epsilon,rc=rc_t_e, wc=wc_t_e))
        elif k[0]=='E' and k[1]=='E':
            table.pair_coeff.set('E', 'E', func=cos_potential, rmin=0,rmax=rc_e_e+wc_e_e,coeff=dict(epsilon=epsilon, rc=rc_e_e,wc=wc_e_e))
        elif k[0]=='H' and k[1]=='N': # add attraction between staple and lipid head
            table.pair_coeff.set('H', 'N', func=cos_potential, rmin=0,rmax=rc_h_n+wc_h_n,coeff=dict(epsilon=0, rc=rc_h_n,wc=wc_h_n))
            #table.pair_coeff.set('H', 'N', func=cos_potential, rmin=0,rmax=rc_h_n+wc_h_n,coeff=dict(epsilon=const_n*epsilon, rc=rc_h_n,wc=wc_h_n))
        #elif k[0]=='H' and k[1]=='N2': # add attraction between staple and lipid head
            #table.pair_coeff.set('H', 'N2', func=cos_potential, rmin=0,rmax=rc_h_n+wc_h_n,coeff=dict(epsilon=const_n*epsilon, rc=rc_h_n,wc=wc_h_n))
        elif k[0]=='H' and k[1]=='C': # add attraction between staple and lipid head
            table.pair_coeff.set('H', 'C', func=cos_potential, rmin=0,rmax=rc_h_c+wc_h_c,coeff=dict(epsilon=0, rc=rc_h_c,wc=wc_h_c))
            #table.pair_coeff.set('H', 'C', func=cos_potential, rmin=0,rmax=rc_h_c+wc_h_c,coeff=dict(epsilon=const_c*epsilon, rc=rc_h_c,wc=wc_h_c))
        else:
            command="table.pair_coeff.set('%s','%s', func=cos_potential,rmin=0,rmax=rc_t_t+wc, coeff=dict(epsilon=0, rc=rc_t_t,wc=wc))"%(k[0],k[1])
            eval(command)

    # Repulsive Shifted LJ Potential
    lj_shift = hoomd.md.pair.lj(10*sigma_bead, nlist=nl) # ignore 10 sigma: it is a large number that is irrelevant
    lj_shift.set_params(mode='shift') # this does the plus epsilon term for lj
    for k in itertools.combinations_with_replacement(p_type_arr,2):
        if k[0] in ['VA','STX','SN1','SN2','LL'] or k[1] in ['VA','STX','SN1','SN2','LL']: 
            command="lj_shift.pair_coeff.set('%s','%s',epsilon=0, sigma=b_%s_%s,r_cut = rc_%s_%s)"%(k[0],k[1],k[0].lower(),k[1].lower(),k[0].lower(),k[1].lower())
            eval(command)
        else:
            if k[0] == 'G' and k[1] == 'G':
                command="lj_shift.pair_coeff.set('%s','%s',epsilon=0, sigma=b_%s_%s,r_cut = rc_%s_%s)"%(k[0],k[1],k[0].lower(),k[1].lower(),k[0].lower(),k[1].lower())
                eval(command)
            elif 'L' in [k[0],k[1]]:
                command="lj_shift.pair_coeff.set('%s','%s',epsilon=0, sigma=b_%s_%s,r_cut = rc_%s_%s)"%(k[0],k[1],k[0].lower(),k[1].lower(),k[0].lower(),k[1].lower())
                eval(command)
            elif 'LJ' in [k[0],k[1]]:
                command="lj_shift.pair_coeff.set('%s','%s',epsilon=0, sigma=b_%s_%s,r_cut = rc_%s_%s)"%(k[0],k[1],k[0].lower(),k[1].lower(),k[0].lower(),k[1].lower())
                eval(command)
            elif 'LI' in [k[0],k[1]]:
                command="lj_shift.pair_coeff.set('%s','%s',epsilon=0, sigma=b_%s_%s,r_cut = rc_%s_%s)"%(k[0],k[1],k[0].lower(),k[1].lower(),k[0].lower(),k[1].lower())
                eval(command)
            elif 'LJI' in [k[0],k[1]]:
                command="lj_shift.pair_coeff.set('%s','%s',epsilon=0, sigma=b_%s_%s,r_cut = rc_%s_%s)"%(k[0],k[1],k[0].lower(),k[1].lower(),k[0].lower(),k[1].lower())
                eval(command)
            elif 'LT' in [k[0],k[1]]:
                command="lj_shift.pair_coeff.set('%s','%s',epsilon=0, sigma=b_%s_%s,r_cut = rc_%s_%s)"%(k[0],k[1],k[0].lower(),k[1].lower(),k[0].lower(),k[1].lower())
                eval(command)  
            else:
                command="lj_shift.pair_coeff.set('%s','%s',epsilon=epsilon, sigma=b_%s_%s,\
                r_cut = rc_%s_%s)" % (k[0],k[1],k[0].lower(),k[1].lower(),k[0].lower(),k[1].lower())
                eval(command)

    # Set SNARE-SNARE electrostatic potential
    # NOT USED - argv[4] flag always set to 0 for SNARE runs
    # r,rmin,rmax,r_bead1,r_bead2,charge1,charge2 
    if int(sys.argv[4])!=0:
        table = hoomd.md.pair.table(width=100000,nlist = nl,name='elect_pot')

        cut_off_electro = int(3)*0.8
    
        for k in itertools.combinations_with_replacement(p_type_arr,2):
            if k[0] in ['n2J','n1J','0J','1J','2J','3J']:
                if k[1] in ['n2J','n1J','0J','1J','2J','3J']:
                    command = "table.pair_coeff.set('%s','%s', func=elect_pot,rmin=0,rmax=cut_off_electro/0.88+0.66/0.88, coeff=dict(r_bead1=b_%s, r_bead2=b_%s, charge1=charge_%s, charge2=charge_%s, wcharge=1))"%(k[0],k[1],k[0].lower(),k[1].lower(),k[0].lower(),k[1].lower())
                    eval(command)
                else:
                    command = "table.pair_coeff.set('%s','%s', func=elect_pot,rmin=0,rmax=cut_off_electro/0.88+0.66/0.88, coeff=dict(r_bead1=b_%s, r_bead2=b_%s, charge1=0, charge2=0, wcharge=0))"%(k[0],k[1],k[0].lower(),k[1].lower())
                    eval(command)
            else:
                command = "table.pair_coeff.set('%s','%s', func=elect_pot,rmin=0,rmax=cut_off_electro/0.88+0.66/0.88, coeff=dict(r_bead1=b_%s, r_bead2=b_%s, charge1=0, charge2=0, wcharge=0))"%(k[0],k[1],k[0].lower(),k[1].lower())
                eval(command)

    # TMD wall potentials (repulsive part of the LJ potential)
    ws = hoomd.md.wall.group()
    epsilon_w = 0.1*epsilon
    sigma_w = 1.5*tmd_bd_r
    delta_r_w = 2*sigma_w + tmd_bd_r/2 # this is the radial distance from wall to wall
    r_cut_w = 0.99*delta_r_w/2
    r_cyl_in = tmd_pos_r_avg - delta_r_w/2
    r_cyl_out = tmd_pos_r_avg + delta_r_w/2

    ws.add_cylinder(r=r_cyl_in, origin=(0,0,0), axis=(0,0,1),inside = False)
    ws.add_cylinder(r=r_cyl_out, origin=(0,0,0), axis=(0,0,1),inside = True)
    wlj = hoomd.md.wall.lj(ws,r_cut = r_cut_w)

    for k in p_type_arr:
        if k=='D' or k == 'D2':
            wlj.force_coeff.set([k],epsilon=epsilon_w,sigma=sigma_w,alpha=0,r_cut=r_cut_w)
        else:
            wlj.force_coeff.set([k],epsilon=0,sigma=sigma_w,alpha=0,r_cut=r_cut_w)

    # Bond Potentials
    fene = hoomd.md.bond.fene()
    fene.bond_coeff.set('bond1', k=k_bond, r0=r_inf, sigma=b_h_t, epsilon=epsilon)
    fene.bond_coeff.set('bond2', k=k_bond, r0=r_inf, sigma=b_t_t, epsilon=epsilon)
    fene.bond_coeff.set('bond3', k=k_bond, r0=r_inf, sigma=b_t_t, epsilon=epsilon)
    fene.bond_coeff.set('bond4', k=0, r0=10*r_inf, sigma=b_h_t, epsilon=epsilon)
    if int(sys.argv[2])==0:
        fene.bond_coeff.set('linker', k=0, r0=10*r_inf, sigma=b_n_n1j, epsilon=0)
        fene.bond_coeff.set('linker2', k=0, r0=10*r_inf, sigma=b_n_n1j, epsilon=0)
    else:
        fene.bond_coeff.set('linker', k=0, r0=10*r_inf, sigma=b_n_l, epsilon=0)
    
    # Bending Potentials
    harmonic = hoomd.md.bond.harmonic()
    harmonic.bond_coeff.set('bond1', k=0, r0=2*sigma_bead) 
    harmonic.bond_coeff.set('bond2', k=0, r0=2*sigma_bead) 
    harmonic.bond_coeff.set('bond3', k=0, r0=2*sigma_bead) 
    harmonic.bond_coeff.set('bond4', k=k_bend, r0=6*sigma_bead)
    harmonic.bond_coeff.set('linker', k=0,r0=0*r_inf)
    harmonic.bond_coeff.set('linker2', k=0,r0=0*r_inf)

    # Linker Potentials
    # sys.argv[3] flag set to 0 for SNARE runs
    # r,rmin,rmax,N_unzip,linker
    global ld_table
    ld_table = hoomd.md.bond.table(width=100000)
    if int(sys.argv[3])==0: 
        ld_table.bond_coeff.set('bond1',func=linker_pot_real_snare,rmin=0,rmax=10,coeff=dict(N_unzip=N_unzip,lp=lp,linker=0))
        ld_table.bond_coeff.set('bond2',func=linker_pot_real_snare,rmin=0,rmax=10,coeff=dict(N_unzip=N_unzip,lp=lp,linker=0))
        ld_table.bond_coeff.set('bond3',func=linker_pot_real_snare,rmin=0,rmax=10,coeff=dict(N_unzip=N_unzip,lp=lp,linker=0))
        ld_table.bond_coeff.set('bond4',func=linker_pot_real_snare,rmin=0,rmax=10,coeff=dict(N_unzip=N_unzip,lp=lp,linker=0))
        ld_table.bond_coeff.set('linker',func=linker_pot_real_snare,rmin=0,rmax=70,coeff=dict(N_unzip=N_unzip,lp=lp,linker=1))
        ld_table.bond_coeff.set('linker2',func=linker_pot_real_snare,rmin=0,rmax=70,coeff=dict(N_unzip=N_unzip,lp=lp,linker=1))
    else:
        ld_table.bond_coeff.set('bond1',func=linker_pot,rmin=0,rmax=12,coeff=dict(k=0,r1=1))
        ld_table.bond_coeff.set('bond2',func=linker_pot,rmin=0,rmax=12,coeff=dict(k=0,r1=1))
        ld_table.bond_coeff.set('bond3',func=linker_pot,rmin=0,rmax=12,coeff=dict(k=0,r1=1))
        ld_table.bond_coeff.set('bond4',func=linker_pot,rmin=0,rmax=12,coeff=dict(k=0,r1=1))
        ld_table.bond_coeff.set('linker',func=linker_pot,rmin=0,rmax=70,coeff=dict(k=0,r1=r1_linker))
    
    ## CREATE GROUPS
    groupH = hoomd.group.type(name='groupH', type='H', update = True)
    groupT = hoomd.group.type(name='groupT', type='T', update = True)
    groupG = hoomd.group.type(name='groupG', type='G', update = True)
    groupROD1 = hoomd.group.type(name='groupROD1', type='A', update = True)
    groupROD2 = hoomd.group.type(name='groupROD2', type='A2', update = True)
    groupROD = hoomd.group.union(name='groupROD', a=groupROD1, b=groupROD2)
    groupTMD1 = hoomd.group.type(name='groupTMD1', type='D', update = True)
    groupTMD2 = hoomd.group.type(name='groupTMD2', type='D2', update = True)
    groupTMD = hoomd.group.union(name='groupTMD', a=groupTMD1, b=groupTMD2)
    groupHT=hoomd.group.union(name='groupHT', a=groupH, b=groupT)
    groupHTG=hoomd.group.union(name='groupHTG',a=groupHT,b=groupG)
    rigid_center=hoomd.group.rigid_center()
    gsd_group=hoomd.group.union(name='gsd_group',a=groupHTG,b=rigid_center)

    global npy_name
    seed = '_%.0f'%seed_sim
    npy_name =f_name+seed
    
    snare_npy_name = f_name+'_snare_force'+seed+'.npy'
    snare_pos_npy_name = f_name+'_snare_pos'+seed+'.npy'
    
    snare_force_arr_all = []
    snare_pos_arr_all = []

    particles_A_ids = [p.tag for p in system.particles if p.type == 'A' or p.type == 'A2']
    particles_A_bodies = [system.particles[pid].body for pid in particles_A_ids]
    print("IDs of SNARE center beads:" ,particles_A_ids)
    print("Bodies of SNARE center beads:" ,particles_A_bodies)

    # Save the moving average force acting on each SNARE COM and its instantaneous position
    forces_A_temp = []
    window = 200
    def log_snare_force(timestep):
        global snare_force_arr_all, snare_pos_arr_all, forces_A_temp
         
        forces_A = [system.particles[pid].net_force for pid in particles_A_ids]
        forces_A_temp.append(forces_A)

        if timestep%window == 0:
            forces_A_avg = np.mean(forces_A_temp, axis=0)
            snare_force_arr_all.append(forces_A_avg)
            forces_A_temp = []

            pos_A = [system.particles[pid].position for pid in particles_A_ids]
            snare_pos_arr_all.append(pos_A)

    # Save the id and force acting on each SNARE linker
    linker_bonds = [b for b in system.bonds if b.type == 'linker' or b.type == 'linker2']

    # Save the ids of all particles involve in linker force
    up_linker_particles_ct = [bond.a for bond in linker_bonds if (system.particles[bond.a].type == 'LJ' or system.particles[bond.a].type == 'LJI') and system.particles[bond.a].position[2] > 0]
    up_linker_particles_ct = up_linker_particles_ct + [bond.b for bond in linker_bonds if (system.particles[bond.b].type == 'LJ' or system.particles[bond.b].type == 'LJI') and system.particles[bond.b].position[2] > 0]
    up_linker_particles_nt = [bond.b for bond in linker_bonds if (system.particles[bond.a].type == 'LJ' or system.particles[bond.a].type == 'LJI') and system.particles[bond.a].position[2] > 0]
    up_linker_particles_nt = up_linker_particles_nt + [bond.a for bond in linker_bonds if (system.particles[bond.b].type == 'LJ' or system.particles[bond.b].type == 'LJI') and system.particles[bond.b].position[2] > 0]
    
    down_linker_particles_ct = [bond.a for bond in linker_bonds if (system.particles[bond.a].type == 'LJ' or system.particles[bond.a].type == 'LJI') and system.particles[bond.a].position[2] < 0]
    down_linker_particles_ct = down_linker_particles_ct + [bond.b for bond in linker_bonds if (system.particles[bond.b].type == 'LJ' or system.particles[bond.b].type == 'LJI') and system.particles[bond.b].position[2] < 0]
    down_linker_particles_nt = [bond.b for bond in linker_bonds if (system.particles[bond.a].type == 'LJ' or system.particles[bond.a].type == 'LJI') and system.particles[bond.a].position[2] < 0]
    down_linker_particles_nt = down_linker_particles_nt + [bond.a for bond in linker_bonds if (system.particles[bond.b].type == 'LJ' or system.particles[bond.b].type == 'LJI') and system.particles[bond.b].position[2] < 0]

    # Iterate over up_linker_particles, and for each pair, print the body attribute of the particle of type '1J'
    up_linker_bodies = []
    for i in up_linker_particles_ct:
        up_linker_bodies.append(system.particles[i].body)
    
    # Iterate over down_linker_particles, and for each pair, print the body attribute of the particle of type '2J'
    down_linker_bodies = []
    for i in down_linker_particles_ct:
        down_linker_bodies.append(system.particles[i].body)

    # Print the bodies of all linker particles
    print("Bodies of up_linker_particles:" ,up_linker_bodies)
    print("Bodies of down_linker_particles:" ,down_linker_bodies)
    '''
    '''
    txt_file_name = f_name+'_body_order.txt'
    with open(txt_file_name, 'w') as file:
        file.write("IDs of SNARE center beads:\n")
        for i in particles_A_ids:
            file.write("%i " % i)
        file.write("\nBodies of SNARE center beads:\n")
        for i in particles_A_bodies:
            file.write("%i " % i)
        file.write("\nBodies of up_linker_particles:\n")
        for i in up_linker_bodies:
            file.write("%i " % i)
        file.write("\nBodies of down_linker_particles:\n")
        for i in down_linker_bodies:
            file.write("%i " % i)
    '''
    '''
    ld_up_force_npy_name =f_name+'_force_up'+seed+'.npy'
    ld_up_dir_npy_name =f_name+'_dir_up'+seed+'.npy'
    ld_up_ct_npy_name =f_name+'_ct_up'+seed+'.npy'
    linker_force_up_all = []
    linker_dir_up_all = []
    linker_ct_up_all = []

    ld_f_up_temp = []
    def log_ld_force_up(timestep):
        global linker_force_up_all, linker_dir_up_all, linker_ct_up_all, ld_f_up_temp
        forces_step = []
        direction_step = []
        ct_step = []
        for i in range(len(up_linker_particles_ct)):
            ct_pos = system.particles[up_linker_particles_ct[i]].position
            nt_pos = system.particles[up_linker_particles_nt[i]].position

            # Calculate the bond length
            dx = ct_pos[0]-nt_pos[0]
            dy = ct_pos[1]-nt_pos[1]
            dz = ct_pos[2]-nt_pos[2]

            r = np.sqrt((dx**2)+(dy**2)+(dz**2)) # in sigma
            V,F = linker_pot_real_snare(r,rmin=0,rmax=70,N_unzip=N_unzip,lp=lp,linker=1)
            F = abs(F)
            forces_x = F*dx/r
            forces_y = F*dy/r
            forces_z = F*dz/r
            forces_step.append([forces_x,forces_y,forces_z])
            if timestep%window == 0:
                direction_step.append((dx,dy,dz))
                ct_step.append(ct_pos)

        ld_f_up_temp.append(forces_step)

        if timestep%window == 0:
            ld_f_up_avg = np.mean(ld_f_up_temp, axis=0)
            linker_force_up_all.append(ld_f_up_avg)
            ld_f_up_temp = []

            linker_dir_up_all.append(direction_step)
            linker_ct_up_all.append(ct_step)

    ld_dn_force_npy_name =f_name+'_force_dn'+seed+'.npy'
    ld_dn_dir_npy_name =f_name+'_dir_dn'+seed+'.npy'
    ld_dn_ct_npy_name =f_name+'_ct_dn'+seed+'.npy'
    linker_force_dn_all = []
    linker_dir_dn_all = []
    linker_ct_dn_all = []

    ld_f_dn_temp = []
    def log_ld_force_dn(timestep):
        global linker_force_dn_all, linker_dir_dn_all, linker_ct_dn_all, ld_f_dn_temp
        forces_step = []
        direction_step = []
        ct_step = []
        for i in range(len(down_linker_particles_ct)):
            ct_pos = system.particles[down_linker_particles_ct[i]].position
            nt_pos = system.particles[down_linker_particles_nt[i]].position

            # Calculate the bond length
            dx = ct_pos[0]-nt_pos[0]
            dy = ct_pos[1]-nt_pos[1]
            dz = ct_pos[2]-nt_pos[2]

            r = np.sqrt((dx**2)+(dy**2)+(dz**2)) # in sigma
            V,F = linker_pot_real_snare(r,rmin=0,rmax=70,N_unzip=N_unzip,lp=lp,linker=1)
            F = abs(F)
            forces_x = F*dx/r
            forces_y = F*dy/r
            forces_z = F*dz/r
            forces_step.append([forces_x,forces_y,forces_z])
            if timestep%window == 0:
                direction_step.append((dx,dy,dz))
                ct_step.append(ct_pos)

        ld_f_dn_temp.append(forces_step)

        if timestep%window == 0:
            ld_f_dn_avg = np.mean(ld_f_dn_temp, axis=0)
            linker_force_dn_all.append(ld_f_dn_avg)
            ld_f_dn_temp = []

            linker_dir_dn_all.append(direction_step)
            linker_ct_dn_all.append(ct_step)
    
    def save_files(timestep):
        global snare_force_arr_all, snare_pos_arr_all, linker_force_up_all, linker_dir_up_all, linker_ct_up_all, linker_force_dn_all, linker_dir_dn_all, linker_ct_dn_all
        np.save(snare_npy_name,snare_force_arr_all)
        np.save(snare_pos_npy_name,snare_pos_arr_all)
        np.save(ld_up_force_npy_name,linker_force_up_all)
        np.save(ld_up_dir_npy_name,linker_dir_up_all)
        np.save(ld_up_ct_npy_name,linker_ct_up_all)
        np.save(ld_dn_force_npy_name,linker_force_dn_all)
        np.save(ld_dn_dir_npy_name,linker_dir_dn_all)
        np.save(ld_dn_ct_npy_name,linker_ct_dn_all)


    ###
    ### Perform Equilibration
    ###
    hoomd.md.integrate.mode_standard(dt=dt_eq)
    bd1 = hoomd.md.integrate.langevin(group=gsd_group,kT=T,seed=seed_sim)
    Gamma_snare = Gamma * 69
    bd1.set_gamma('H',Gamma)
    bd1.set_gamma('T',Gamma)
    bd1.set_gamma('G',GammaG)
    bd1.set_gamma('D',Gamma)
    bd1.set_gamma('D2',Gamma)
    bd1.set_gamma('A',Gamma_snare)
    bd1.set_gamma('A2',Gamma_snare)
    hoomd.compute.thermo(group=groupROD)
    N_eq=1500000 # ~1 microsec equilibration
    seed = '_%.0f'%seed_sim
    #gsd_name_equ = f_name+seed+'_equ'+'.gsd'
    #movie_equ=hoomd.dump.gsd(gsd_name_equ,period=20000,group=hoomd.group.all(),overwrite=True)
    hoomd.run(N_eq)

    

    ###
    ### Run simulation
    ###
    #movie_equ.disable()
    force_period = 10
    ld_force_up_logger = hoomd.analyze.callback(callback=log_ld_force_up,period=force_period)
    ld_force_dn_logger = hoomd.analyze.callback(callback=log_ld_force_dn,period=force_period)
    snare_force_logger = hoomd.analyze.callback(callback=log_snare_force, period=force_period)
    hoomd.md.integrate.mode_standard(dt=dt)
    bd1.set_gamma('H',Gamma)
    bd1.set_gamma('T',Gamma)
    bd1.set_gamma('G',GammaG)
    bd1.set_gamma('D',Gamma)
    bd1.set_gamma('D2',Gamma)
    bd1.set_gamma('A',Gamma_snare)
    bd1.set_gamma('A2',Gamma_snare)
    N_run = 30000000 # number of timesteps, 1 timestep = 0.068 ns
    gsd_period = 100000
    gsd_name = f_name+seed+'.gsd'
    hoomd.dump.gsd(gsd_name,period=gsd_period,group=hoomd.group.all(),overwrite=True)
    force_log_period = 1000000
    save_files_logger = hoomd.analyze.callback(callback=save_files, period=force_log_period)
    gsd_check_name = f_name+seed+'_chk'+'.gsd'
    hoomd.run(N_run)

