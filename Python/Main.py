"""
# Main code for Sufficient Statistic project
"""
import sys 
sys.path.insert(0,'OneAssetBond/')

import numpy as np
import defineSSParametersIOUsBond as Params
from copy import copy
import pickle
from SteadyState_fsolve import SteadyState_fsolve
from FluctuationsOneAssetIOUsBond import FluctuationsOneAssetIOUsBond, Fsys, SGU_solver, plot_IRF

EconomyParams = copy(Params.parm_one_asset_IOUsBond)
  

SSEconomy = SteadyState_fsolve(**EconomyParams)


##### Choose whether calculate new steady state or use ond one

#SSS = SSEconomy.SolveSteadyState() # New steady state
#pickle.dump(SSS, open("btoy5_tau_1.p", "wb"))

SSS=pickle.load(open("btoy5_tau_1.p", "rb")) # Use old steady state

#############3#################################################3



##### Preparation: Calculation of sufficient statistics

#   MPCs
par = SSS['par']
mpar = SSS['mpar']
targets = SSS['targets']
grid = SSS['grid']
c_guess = SSS['c_policy']
inc = SSS['inc']
jd = SSS['joint_distr'].copy()
joint_distr = jd.copy()
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
#    MPC_m[:,hh]=np.gradient(np.squeeze(C_ind[:,hh]))/np.gradient(grid['m'])  # MPC_m_ is same with MPC_m
    MPC_m[:,hh]=np.gradient(np.squeeze(c_guess[:,hh]))/np.gradient(grid['m'])  # MPC_m_ is same with MPC_m
MPC_h = np.zeros((mpar['nm'],mpar['nh']))

for mm in range(0,mpar['nm']):
   MPC_h[mm,:]=np.gradient(C_ind[mm,:])/np.gradient(np.multiply(WW_h,grid_h_aux))

NNP = meshes_m
URE = inc['labor'] + meshes_m - c_guess + inc['profits']

# MPC_m=np.min(MPC_m,1)

##### Calculation of sufficient statistics

# 1. Aggregate income channel: EI[Yi/Y MPC_m]
Inc_wt_MPC = np.sum(np.sum( np.multiply(np.multiply(WW_h_mesh/np.sum(np.sum(np.multiply(jd.copy(),WW_h_mesh))), MPC_m), jd.copy()) ))

# 3. Fisher channel
Redist_elas_P = np.sum(np.sum(np.multiply(np.multiply(MPC_m.copy(), NNP.copy()), jd.copy()))) - np.sum(np.sum(np.multiply(MPC_m.copy(),jd.copy())))*np.sum(np.sum(np.multiply(NNP.copy(),jd.copy())))

# 4. Interest rate exposure channel
Redist_elas_R = np.sum(np.sum(np.multiply(np.multiply(MPC_m.copy(),URE.copy()),jd.copy()))) - np.sum(np.sum(np.multiply(MPC_m.copy(),jd.copy())))*np.sum(np.sum(np.multiply(URE.copy(),jd.copy())))

# 5. substitution channel
#sig_i = par['xi']**(-1)*np.multiply(c_guess,1/C_ind)
#Hick_scaling = np.sum(np.sum( np.multiply(np.multiply(np.multiply(sig_i, (1-MPC_m.copy())), C_ind.copy()), jd.copy()) ))

sig_i = par['xi']**(-1)
Hick_scaling = np.sum(np.sum( np.multiply(np.multiply(np.multiply(sig_i, (1-MPC_m.copy())), c_guess.copy()), jd.copy()) ))


# 6. additional term for GHH preference



#### changes of policy parameters: irrelevant of a steady state

#SSS['par']['rho_B'] = 0.99
#SSS['par']['gamma_pi'] = 1.25
SSS['par']['theta_pi'] = 2
##############################################################################

#pickle.dump(SSS, open("SSS.p", "wb"))
#
#SSS=pickle.load(open("SSS.p", "rb"))


EX2SR=FluctuationsOneAssetIOUsBond(**SSS)

SR=EX2SR.StateReduc()

State       = np.zeros((SR['mpar']['numstates'],1))
State_m     = State
Contr       = np.zeros((SR['mpar']['numcontrols'],1))
Contr_m     = Contr

Fsysresult = Fsys(State, State_m, Contr, Contr_m, SR['Xss'], SR['Yss'], 
         SR['Gamma_state'], SR['Gamma_control'], SR['InvGamma'], SR['Copula'],SR['par'], 
         SR['mpar'], SR['grid'], SR['targets'],SR['P_H'],SR['aggrshock'],SR['oc'])
    
    
SGUresult=SGU_solver(SR['Xss'],SR['Yss'],SR['Gamma_state'],SR['Gamma_control'],SR['InvGamma'],SR['Copula'],
                     SR['par'],SR['mpar'],SR['grid'],SR['targets'],SR['P_H'],SR['aggrshock'],SR['oc'])


plotresult = plot_IRF(SR['mpar'],SR['par'],SGUresult['gx'],SGUresult['hx'],SR['joint_distr'],
         SR['Gamma_state'],SR['grid'],SR['targets'],SR['os'],SR['oc'],SR['Output'],C_ind,c_guess,
         WW_h_mesh,MPC_m,Inc_wt_MPC,Redist_elas_P,Redist_elas_R,Hick_scaling,URE,NNP)

IRF_W = plotresult['IRF_W']
IRF_N = plotresult['IRF_N']
IRF_Profit = plotresult['IRF_Profit']
    
print(round(plotresult['IRF_Xagg'][0,0],3))
print(round(plotresult['IRF_X_by_suff'][0,0],3))
