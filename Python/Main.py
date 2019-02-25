"""
# Main code for Sufficient Statistic project
"""
import sys 
sys.path.insert(0,'OneAssetBond/')

import numpy as np
import defineSSParameters as Params
from copy import copy
import pickle
from SteadyStateOneAssetIOUsBond import SteadyStateOneAssetIOUsBond


EconomyParams = copy(Params.parm_one_asset_IOUsBond)
  
SSEconomy = SteadyStateOneAssetIOUsBond(**EconomyParams)

SSSolution = SSEconomy.SolveSteadyState()

pickle.dump(SSSolution, open("SSSolution.p", "wb"))

SSSolution=pickle.load(open("SSSolution.p", "rb"))

#   MPCs
par = SSSolution['par']
mpar = SSSolution['mpar']

##  Calculate MPC
MPC_m = np.zeros((mpar['nm'],mpar['nh']))
for hh in range(mpar['nh']):
    MPC_m[:,hh]=np.gradient(SSSolution['c_policy_noncomposite'][:,hh])/np.gradient(SSSolution['grid']['m'])
    
NNP = SSSolution['meshesm']
# Note: 'labor' income is for the composite, but URE excludes leisure spending. 
# That is why we take away 'c_policy' which is spending on the composite
# good, not 'c_policy_noncomposite'
URE = SSSolution['inc']['labor'] + SSSolution['meshesm'] - SSSolution['c_policy']

## Calculate Sufficient statistics in Auclert
#
agg_c = sum(sum(SSSolution['joint_distr']*SSSolution['c_policy_noncomposite']))
# 1. Aggregate income channel: EI[Yi/Y MPC_m]
# income includes interest rate income?
Inc_wt_MPC = sum(sum( SSSolution['labor_income_noncomposite']*MPC_m*SSSolution['joint_distr'] ))/agg_c
#
# 3. Fisher channel
Redist_elas_P = sum(sum(MPC_m*NNP*SSSolution['joint_distr']))/agg_c - sum(sum(MPC_m*SSSolution['joint_distr']))*sum(sum(NNP*SSSolution['joint_distr']))/agg_c
#
# 4. Interest rate exposure channel
Redist_elas_R = sum(sum(MPC_m*URE*SSSolution['joint_distr']))/agg_c - sum(sum(MPC_m*SSSolution['joint_distr']))*sum(sum(URE*SSSolution['joint_distr']))/agg_c
#
# 5. substitution channel - note in our GHH model the elasticity of 
# intertemporal substitution is different for different households
sig_i = par['xi']**(-1)*SSSolution['c_policy']/SSSolution['c_policy_noncomposite'] # c_policy: composite goods consumption, c_policy_noncomposite: final goods consumption
#
Intertemporal_subs = sum(sum(sig_i*(1.0-MPC_m)*SSSolution['c_policy_noncomposite']*SSSolution['joint_distr']))/agg_c

