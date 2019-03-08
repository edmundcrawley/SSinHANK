"""
# Main code for Sufficient Statistic project
"""
import sys 
sys.path.insert(0,'OneAssetBond/')

import numpy as np
import defineSSParametersIOUsBond as Params
from copy import copy
import pickle
from SteadyStateOneAssetIOUsBond import SteadyStateOneAssetIOUsBond
from FluctuationsOneAssetIOUsBond import FluctuationsOneAssetIOUsBond, SGU_solver, plot_IRF

EconomyParams = copy(Params.parm_one_asset_IOUsBond)

        
SSEconomy = SteadyStateOneAssetIOUsBond(**EconomyParams)

SSS = SSEconomy.SolveSteadyState()


##### Preparation: Calculation of sufficient statistics

#   MPCs
par = SSS['par']
mpar = SSS['mpar']
targets = SSS['targets']
grid = SSS['grid']
c_guess = SSS['c_policy']
inc = SSS['inc']
jd = SSS['joint_distr'].copy()
NW=(targets['N']/par['H'])*targets['W']
WW=NW*np.ones((mpar['nm'],mpar['nh'])) #Wages
WW[:,-1]=SSS['Profits_fc']*par['profitshare']

meshes_m,meshes_h = np.meshgrid(grid['m'],grid['h'], indexing='ij')
aux_x = par['tau']*targets['N']/par['H']*targets['W']*meshes_h/(1+par['gamma'])
aux_x[:,-1]=0
C_ind = c_guess + aux_x


# MPC
WW_h=WW[0,:]
WW_h_mesh=np.multiply(WW,meshes_h)
grid_h_aux=grid['h']
MPC_m = np.zeros((mpar['nm'],mpar['nh']))

for hh in range(0,mpar['nh']):
    MPC_m[:,hh]=np.gradient(np.squeeze(C_ind[:,hh]))/np.gradient(grid['m'])  # MPC_m_ is same with MPC_m

MPC_h = np.zeros((mpar['nm'],mpar['nh']))

for mm in range(0,mpar['nm']):
   MPC_h[mm,:]=np.gradient(C_ind[mm,:])/np.gradient(np.multiply(WW_h,grid_h_aux))

NNP = meshes_m
URE = inc['labor'] + meshes_m - c_guess + inc['profits']

# MPC_m=np.min(MPC_m,1)

##### Calculation of sufficient statistics

# 1. Aggregate income channel: EI[Yi/Y MPC_m]
Inc_wt_MPC = np.sum(np.sum( np.multiply(WW_h_mesh/np.sum(np.sum(np.multiply(jd.copy(),WW_h_mesh))), MPC_m, jd.copy()) ))

# 3. Fisher channel
Redist_elas_P = np.sum(np.sum(np.multiply(MPC_m.copy(), NNP.copy(), jd.copy()))) - np.sum(np.sum(np.multiply(MPC_m.copy(),jd.copy())))*np.sum(np.sum(np.multiply(NNP,jd.copy())))

# 4. Interest rate exposure channel
Redist_elas_R = np.sum(np.sum(np.multiply(MPC_m.copy(),URE.copy(),jd.copy()))) - np.sum(np.sum(np.multiply(MPC_m.copy(),jd.copy())))*np.sum(np.sum(np.multiply(URE.copy(),jd.copy())))

# 5. substitution channel
sig_i = par['xi']**(-1)*np.multiply(c_guess,1/C_ind)

Hick_scaling = np.sum(np.sum( np.multiply(np.multiply(sig_i, (1-MPC_m.copy()), C_ind.copy()), jd.copy()) ))


#### changes of policy parameters: irrelevant of a steady state

SSS['par']['rho_B'] = 0.99
SSS['par']['gamma_pi'] = 1.25

####

pickle.dump(SSS, open("SSS.p", "wb"))

SSS=pickle.load(open("SSS.p", "rb"))


EX2SR=FluctuationsOneAssetIOUsBond(**SSS)

SR=EX2SR.StateReduc()
    
SGUresult=SGU_solver(SR['Xss'],SR['Yss'],SR['Gamma_state'],SR['Gamma_control'],SR['InvGamma'],SR['Copula'],
                     SR['par'],SR['mpar'],SR['grid'],SR['targets'],SR['P_H'],SR['aggrshock'],SR['oc'])

plot_IRF(SR['mpar'],SR['par'],SGUresult['gx'],SGUresult['hx'],SR['joint_distr'],
         SR['Gamma_state'],SR['grid'],SR['targets'],SR['os'],SR['oc'],SR['Output'])
    
