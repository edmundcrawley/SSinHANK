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
#
##   MPCs
#par = SSSolution['par']
#mpar = SSSolution['mpar']
#targets = SSSolution['targets']
#NW=(targets['N']/par['H'])*targets['W']
#WW=NW*np.ones((mpar['nm'],mpar['nh'])) #Wages
#WW(:,end)=par.PROFITS*par.profitshare;
#
#% MPC
#WW_h=squeeze(WW(1,:)); WW_h_mesh=squeeze(WW(:,:).*meshes.h);
#grid_h_aux=grid.h;
#MPC_m = zeros(mpar.nm,mpar.nh);
#
#for hh=1:mpar.nh
#    MPC_m(:,hh)=gradient(squeeze(C_ind(:,hh)))./gradient(grid.m)'; % MPC_m_ is same with MPC_m
#end
