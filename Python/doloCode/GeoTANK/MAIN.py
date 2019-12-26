import numpy as np
from matplotlib import pyplot as plt
from dolo import *  # load the dolo library
from pathlib import Path

figure_dir = Path('./Figures/')
table_dir = Path('./Tables/')

#################################################
# Geo-RANK model
#################################################          
GeoRANKmodel = yaml_import("GeoRANK.yaml")   # import the model
T_sim = 20
GeoRANK_pert = approximate_controls(GeoRANKmodel, order=1)
GeoRANK_irf = response(GeoRANKmodel,GeoRANK_pert, 'eps_G', T_sim)

plt.plot(GeoRANK_irf.sel(V='dC'))
plt.plot(GeoRANK_irf.sel(V='dY'))
plt.plot(GeoRANK_irf.sel(V='dG'))
#plt.plot(GeoRANK_irf.sel(V='dT'))
#plt.plot(GeoRANK_irf.sel(V='dA'))
#plt.plot(irf.sel(V='dA'))
plt.legend(loc='upper right', labels = ["dC","dY","dG"])
plt.savefig(figure_dir / 'GeoRANK_irf.pdf')
plt.show()
#################################################
# Geo-TANK model
#################################################  

GeoTANKmodel = yaml_import("GeoTANK.yaml")   # import the model

GeoTANK_pert = approximate_controls(GeoTANKmodel, order=1)
GeoTANK_irf = response(GeoTANKmodel,GeoTANK_pert, 'eps_G', T_sim)

plt.plot(GeoTANK_irf.sel(V='dC'))
plt.plot(GeoTANK_irf.sel(V='dY'))
plt.plot(GeoTANK_irf.sel(V='dG'))
#plt.plot(GeoTANK_irf.sel(V='dT'))
#plt.plot(GeoTANK_irf.sel(V='dA'))
plt.legend(loc='upper right', labels = ["dC","dY","dG"])
plt.savefig(figure_dir / 'GeoTANK_irf.pdf')
plt.show()

length= T_sim
plt.plot(GeoTANK_irf.sel(V='dC')[0:length])
plt.plot(GeoRANK_irf.sel(V='dC')[0:length])
plt.show()

#################################################
num_points = 100
T_sim = 4000
this_rho_B = np.linspace(0,0.9999999,num_points)

impact_multiplier_TANK = np.zeros(num_points)
cum_multiplier_TANK = np.zeros(num_points)
for i in range(num_points):
    GeoTANKmodel.set_calibration(rho_B=this_rho_B[i])
    GeoTANK_pert = approximate_controls(GeoTANKmodel, order=1)
    GeoTANK_irf = response(GeoTANKmodel,GeoTANK_pert, 'eps_G', 4000)
    impact_multiplier_TANK[i] = GeoTANK_irf.sel(V='dY')[2]
    cum_Y_TANK = np.sum(1.0/1.01**np.arange(T_sim)*GeoTANK_irf.sel(V='dY'))
    cum_G_TANK = np.sum(1.0/1.01**np.arange(T_sim)*GeoTANK_irf.sel(V='dG'))
    cum_multiplier_TANK[i] = cum_Y_TANK/cum_G_TANK
    
impact_multiplier_RANK = np.zeros(num_points)
cum_multiplier_RANK = np.zeros(num_points)
for i in range(num_points):
    GeoRANKmodel.set_calibration(rho_B=this_rho_B[i])
    GeoRANK_pert = approximate_controls(GeoRANKmodel, order=1)
    GeoRANK_irf = response(GeoRANKmodel,GeoRANK_pert, 'eps_G', 4000)
    impact_multiplier_RANK[i] = GeoRANK_irf.sel(V='dY')[2] 
    cum_Y_RANK = np.sum(1.0/1.01**np.arange(T_sim)*GeoRANK_irf.sel(V='dY'))
    cum_G_RANK = np.sum(1.0/1.01**np.arange(T_sim)*GeoRANK_irf.sel(V='dG'))
    cum_multiplier_RANK[i] = cum_Y_RANK/cum_G_RANK
    
T_sim = 4000
GeoTANKmodel.set_calibration(rho_B=0.9)
GeoTANK_pert = approximate_controls(GeoTANKmodel, order=1)
GeoTANK_irf = response(GeoTANKmodel,GeoTANK_pert, 'eps_G', T_sim)
cum_Y_TANK = np.sum(1.0/1.01**np.arange(T_sim)*GeoTANK_irf.sel(V='dY'))
cum_G_TANK = np.sum(1.0/1.01**np.arange(T_sim)*GeoTANK_irf.sel(V='dG'))
cum_multiplier_TANK = cum_Y_TANK/cum_G_TANK
cum_C_TANK = np.sum(1.0/1.01**np.arange(T_sim)*GeoTANK_irf.sel(V='dC'))


GeoRANKmodel.set_calibration(rho_B=0.9)
GeoRANK_pert = approximate_controls(GeoRANKmodel, order=1)
GeoRANK_irf = response(GeoRANKmodel,GeoRANK_pert, 'eps_G', T_sim)
cum_Y_RANK = np.sum(1.0/1.01**np.arange(T_sim)*GeoRANK_irf.sel(V='dY'))
cum_G_RANK = np.sum(1.0/1.01**np.arange(T_sim)*GeoRANK_irf.sel(V='dG'))
cum_multiplier_RANK = cum_Y_RANK/cum_G_RANK
cum_C_RANK = np.sum(1.0/1.01**np.arange(T_sim)*GeoRANK_irf.sel(V='dC'))


plt.plot(this_rho_B,impact_multiplier_RANK)
plt.plot(this_rho_B,impact_multiplier_TANK)
plt.show()

plt.plot(this_rho_B,cum_multiplier_RANK)
plt.plot(this_rho_B,cum_multiplier_TANK)
plt.show()


###########################################################################
# Plot iMPCs
impactMPC = 0.25
T_sim = 30
iMPC = np.zeros((T_sim,T_sim))
for i in range(T_sim):
    iMPC[0,i] = impactMPC*(1-impactMPC)**i
    for j in range(1,T_sim-i):
        iMPC[j,i+j] = iMPC[0,i] 
plt.plot(iMPC[0,:])
plt.plot(iMPC[5,:])
plt.plot(iMPC[10,:])
plt.plot(iMPC[15,:])
plt.savefig(figure_dir / 'iMPC_RANK.pdf')
plt.show()

