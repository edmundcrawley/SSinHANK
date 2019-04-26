import numpy as np
from matplotlib import pyplot as plt
from dolo import *  # load the dolo library
from CalcTransChannels import CalcTransChannels_TANK, CalcTransChannels_TANKcapital
figure_dir = '.\\Figures\\'
table_dir = '.\\Tables\\'

#################################################
# First do the TANK model without capital
#################################################          
TANKmodel = yaml_import("TANKmodel.yaml")   # import the model
#display(TANKmodel)                               # display the model
#TANKmodel.residuals()

param_to_vary = ['pop_share_K','debt_limit']
plot_xlabel = ['Proportion of Keynesian Households','Keynesian Household Debt Level']
param_loop = dict()
param_loop[param_to_vary[0]] = np.linspace(0.0,0.3,31)
param_loop[param_to_vary[1]] = np.linspace(0.0,2.0,41) 

# Loop over sigma
for this_sigma in [1.0,2.0,3.0]:
    TANKmodel.set_calibration(sigma=this_sigma)
    # Loop over the parameter that we are varying
    for j in range(len(param_to_vary)):
        # Reset all varied parameters to their baseline values
        TANKmodel.set_calibration(pop_share_K=0.2)
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
            Transmission_Channels[i,:], suff_stats[i,:], YRP_changes[i,:], checks[i,:] = CalcTransChannels_TANK(TANKmodel)
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
            plt.axvline(x=0.2, color='k', linestyle='--')
            save_name = 'ProportionKeynesian_sigma'+ '{:.0f}'.format(this_sigma) + '.pdf'
        elif param_to_vary[j]=='debt_limit':
            plt.title('$\sigma$ = ' + '{:.1f}'.format(this_sigma))
            save_name = 'KeynesianDebt_sigma'+ '{:.0f}'.format(this_sigma) + '.pdf'
        else:
            print('Can\'t handle this parameter')
        plt.savefig(figure_dir + save_name)


#################################################
# Next do the TANK model WITH capital
#################################################          
TANKcapitalmodel = yaml_import("TANKcapitalmodel.yaml")   # import the model
#display(TANKcapitalmodel)                               # display the model
#TANKcapitalmodel.residuals()

psi_loop = [0.,1.,3.,10000.]
Transmission_Channels = np.zeros((len(psi_loop),6))
suff_stats = np.zeros((len(psi_loop),4))
YRP_changes = np.zeros((len(psi_loop),3))
IRF_i = np.zeros((len(psi_loop),6))
IRF_c_K = np.zeros((len(psi_loop),6))
IRF_c_R = np.zeros((len(psi_loop),6))
IRF_k = np.zeros((len(psi_loop),6))
IRF_q = np.zeros((len(psi_loop),6))
IRF_r_real = np.zeros((len(psi_loop),6))
checks = np.zeros((len(psi_loop),3))

for i in range(len(psi_loop)):
    TANKcapitalmodel.set_calibration(psi_c=psi_loop[i])
    Transmission_Channels[i,:], suff_stats[i,:], YRP_changes[i,:], checks[i,:], \
    IRF_i[i,:], IRF_c_R[i,:], IRF_c_K[i,:], IRF_k[i,:], IRF_q[i,:], IRF_r_real[i,:] \
    = CalcTransChannels_TANKcapital(TANKcapitalmodel)

#Plot IRF for nominal interest rate
plt.figure()
plt.plot(np.transpose(IRF_i), linewidth=1.5)
leg = plt.legend(['$\psi_c = 0$','$\psi_c = 1$','$\psi_c = 3$','$\psi_c = \infty$'], loc='lower right')
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Time (years)')
plt.title('Nominal Interest Rate Path')
plt.ylabel('Nominal Interest Rate')
plt.savefig(figure_dir + 'TANK_capital_IRF_i.pdf')

#Plot IRF for c_K
plt.figure()
plt.plot(np.transpose(IRF_c_K), linewidth=1.5)
leg = plt.legend(['$\psi_c = 0$','$\psi_c = 1$','$\psi_c = 3$','$\psi_c = \infty$'], loc='upper right')
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Time (years)')
plt.title('Keynesian Consumption Path')
plt.ylabel('Consumption')
plt.savefig(figure_dir + 'TANK_capital_IRF_c_K.pdf')

#Plot IRF for c_R
plt.figure()
plt.plot(np.transpose(IRF_c_R), linewidth=1.5)
leg = plt.legend(['$\psi_c = 0$','$\psi_c = 1$','$\psi_c = 3$','$\psi_c = \infty$'], loc='upper right')
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Time (years)')
plt.title('Ricardian Consumption Path')
plt.ylabel('Consumption')
plt.savefig(figure_dir + 'TANK_capital_IRF_c_R.pdf')

#Plot IRF for capital
plt.figure()
plt.plot(np.transpose(IRF_k), linewidth=1.5)
leg = plt.legend(['$\psi_c = 0$','$\psi_c = 1$','$\psi_c = 3$','$\psi_c = \infty$'], loc='upper right')
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Time (years)')
plt.title('Capital Path')
plt.ylabel('Capital')
plt.savefig(figure_dir + 'TANK_capital_IRF_k.pdf')

#Plot IRF for Tobin's q
plt.figure()
plt.plot(np.transpose(IRF_q), linewidth=1.5)
leg = plt.legend(['$\psi_c = 0$','$\psi_c = 1$','$\psi_c = 3$','$\psi_c = \infty$'], loc='upper right')
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Time (years)')
plt.title('Tobin\'s q Path')
plt.ylabel('Tobin\'s q')
plt.savefig(figure_dir + 'TANK_capital_IRF_q.pdf')

#Plot IRF for real interest rate
plt.figure()
plt.plot(np.transpose(IRF_r_real), linewidth=1.5)
leg = plt.legend(['$\psi_c = 0$','$\psi_c = 1$','$\psi_c = 3$','$\psi_c = \infty$'], loc='lower right')
leg.get_frame().set_linewidth(0.0)
plt.xlabel('Time (years)')
plt.title('Real Interest Rate Path')
plt.ylabel('Real Interest Rate')
plt.savefig(figure_dir + 'TANK_capital_IRF_r_real.pdf')

# Now make a table of errors
error_table = '  \\begin{table}\n'
error_table += '\\begin{center}\n'
error_table += '    \\caption{Percentage Error of Decomposition}\\label{table:error}\n'
error_table += '\\begin{tabular}{cccc}  \n'
error_table += '$\\psi_c$ & Total Consumption & Ricardian Consumption & Keynesian Consumption \n'
error_table += '\\\\ \\toprule  \n'
error_table += '{:.0f}'.format(psi_loop[0])+' & '+ '{:.1f}'.format(checks[0,0]) +' \\%  & '+ '{:.1f}'.format(checks[0,1])  +' \\% & 0.0 \\% \\\\ \n'
error_table += '{:.0f}'.format(psi_loop[1])+' & '+ '{:.1f}'.format(checks[1,0]) +' \\%  & '+ '{:.1f}'.format(checks[1,1])  +' \\% & 0.0 \\% \\\\ \n'
error_table += '{:.0f}'.format(psi_loop[2])+' & '+ '{:.1f}'.format(checks[2,0]) +' \\%  & '+ '{:.1f}'.format(checks[2,1])  +' \\% & 0.0 \\% \\\\ \n'
error_table += '$\\infty$ &  '+ '{:.1f}'.format(checks[3,0]) +' \\%  & '+ '{:.1f}'.format(checks[3,1])  +' \\% & 0.0 \\% \\\\ \n'
error_table += '\\\\ \\bottomrule \n '
error_table += '\\end{tabular}\n'
error_table += '\\end{center}\n'
error_table += '\\end{table}\n'
with open(table_dir + 'error_table.tex','w') as f:
    f.write(error_table)
    f.close()
    
# Also write calibration table
param_names = TANKcapitalmodel.symbols['parameters']
param_vals = TANKcapitalmodel.calibration['parameters']
calib_table = '  \\begin{table}\n'
calib_table += '\\begin{center}\n'
calib_table += '    \\caption{Baseline Calibration}\\label{table:calibration}\n'
calib_table += '\\begin{tabular}{lcl}  \n'
calib_table += '\\\\ \\toprule  \n'
calib_table += '$\\sigma$ & '+ '{:.1f}'.format(param_vals[param_names.index('sigma')]) +' & Inverse EIS \\\\ \n'
calib_table += '$\\psi$ & '+ '{:.1f}'.format(param_vals[param_names.index('phi')]) +' & Inverse Frisch Elasticity \\\\ \n'
calib_table += '$\\phi_{\\pi}$ & '+ '{:.1f}'.format(param_vals[param_names.index('phi_pi')]) +' & Taylor Rule Coefficient \\\\ \n'
calib_table += '$\\theta$ & '+ '{:.3f}'.format(param_vals[param_names.index('theta')]) +' & Calvo stickiness parameter \\\\ \n'
calib_table += '$\\beta$ & '+ '{:.1f}'.format(param_vals[param_names.index('beta')]) +' &  Discount Factor\\\\ \n'
calib_table += '$\\alpha$ & '+ '{:.2f}'.format(param_vals[param_names.index('alpha')]) +' &  Capital Share \\\\ \n'
calib_table += '$\\varepsilon$ & '+ '{:.1f}'.format(param_vals[param_names.index('epsilon')]) +' &  Elasticity of sub. between goods\\\\ \n'
calib_table += '$\\lambda$ & '+ '{:.1f}'.format(param_vals[param_names.index('pop_share_K')]) +' & Share of Keynesian Households \\\\ \n'
calib_table += '$\\Omega$ & '+ '{:.1f}'.format(param_vals[param_names.index('debt_limit')]) +' & Keynesian Debt as Share of Income \\\\ \n'
calib_table += '$\\delta$ & '+ '{:.1f}'.format(param_vals[param_names.index('delta')]) +' & Depreciation (capital model only) \\\\ \n'
calib_table += '\\\\ \\bottomrule \n '
calib_table += '\\end{tabular}\n'
calib_table += '\\end{center}\n'
calib_table += '\\end{table}\n'
with open(table_dir + 'calib_table.tex','w') as f:
    f.write(calib_table)
    f.close()


