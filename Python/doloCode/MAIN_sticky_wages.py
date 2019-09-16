import numpy as np
from matplotlib import pyplot as plt
from dolo import *  # load the dolo library
from CalcTransChannels import CalcTransChannels_persistent
from pathlib import Path
figure_dir = Path('./Figures/')
table_dir = Path('./Tables/')

#################################################
# First do the TANK model without capital
#################################################          
TANKmodel = yaml_import("TANK_sticky_wages.yaml")   # import the model
#display(TANKmodel)                               # display the model
#TANKmodel.residuals()

param_to_vary = ['pop_share_K','debt_limit']
plot_xlabel = ['Proportion of Keynesian Households','Keynesian Household Debt Level']
param_loop = dict()
param_loop[param_to_vary[0]] = np.linspace(0.0,0.9,31)
param_loop[param_to_vary[1]] = np.linspace(0.0,2.0,41) 

# Loop over sigma
for this_sigma in [1.0,2.0,3.0,10.0]:
    TANKmodel.set_calibration(sigma=this_sigma)
    # Loop over the parameter that we are varying
    for j in range(len(param_to_vary)):
        # Reset all varied parameters to their baseline values
        TANKmodel.set_calibration(pop_share_K=0.5)
        TANKmodel.set_calibration(debt_limit=0.0)
        param_array = param_loop[param_to_vary[j]]
        Transmission_Channels = np.zeros((len(param_array),6))
        suff_stats = np.zeros((len(param_array),4))
        YRP_changes = np.zeros((len(param_array),3))
        checks = np.zeros((len(param_array),3))
        for i in range(len(param_array)):
            # Set the parameter value - there may be a neater way to do this
            if param_to_vary[j]=='pop_share_K':
                TANKmodel.set_calibration(pop_share_K=param_array[i])
            elif param_to_vary[j]=='debt_limit':
                TANKmodel.set_calibration(debt_limit=param_array[i])
            else:
                print('Can\'t handle this parameter')
            Transmission_Channels[i,:], suff_stats[i,:], YRP_changes[i,:], checks[i,:], IRF_i, IRF_c_R, IRF_c_K, IRF_pi_p, IRF_r_real, IRF_pi_w = CalcTransChannels_persistent(TANKmodel, sticky_wages=True)
        plt.figure()
        plt.stackplot(param_array,np.transpose(Transmission_Channels[:,[4,2,3,0,1]]), \
                      labels=['Intertemporal Substitution', 'Interest Rate Exposure', 'Fisher','Aggregate Income', 'Heterogenous Income'], \
                      edgecolor='k')
        plt.xlabel(plot_xlabel[j])
        plt.ylim([0,3])
        plt.xlim([param_array[0], param_array[-1]])
        plt.ylabel('Consumption Change')
        leg = plt.legend(loc='upper left')
        leg.get_frame().set_linewidth(0.0)
        if param_to_vary[j]=='pop_share_K':
            plt.title('Transmission Channels to a 1% Interest Rate Decline')
            plt.axvline(x=0.5, color='k', linestyle='--')
            save_name = 'ProportionKeynesian_sigma'+ '{:.0f}'.format(this_sigma) + '_sw.pdf'
        elif param_to_vary[j]=='debt_limit':
            plt.title('$\sigma$ = ' + '{:.1f}'.format(this_sigma))
            save_name = 'KeynesianDebt_sigma'+ '{:.0f}'.format(this_sigma) + '_sw.pdf'
        else:
            print('Can\'t handle this parameter')
        plt.savefig(figure_dir + save_name)
