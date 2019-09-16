# Plots delayed response of the average interest rate paid in a given scenario

import numpy as np
from matplotlib import pyplot as plt

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n-1:] / n

bond_quarters = 8
total_quarters = 30
points_per_quarter = 1000
time = np.linspace(0,total_quarters,points_per_quarter*total_quarters)
int_rate_prev = 4.0*np.ones_like(time)
#int_rate_new = int_rate_prev - np.ones_like(time)*np.exp(-time/3.0)

int_rate_new = np.copy(int_rate_prev)
int_rate_new[0:12*points_per_quarter] = 3.0

long_term_rate = moving_average(int_rate_new,bond_quarters*points_per_quarter)

avg_rate = np.zeros_like(long_term_rate)
for i in range(bond_quarters*points_per_quarter):
    avg_rate[i] = (np.sum(long_term_rate[0:i]) + np.sum(int_rate_prev[i+1:bond_quarters*points_per_quarter]))/(bond_quarters*points_per_quarter)
avg_rate[bond_quarters*points_per_quarter:] = moving_average(long_term_rate[1:],bond_quarters*points_per_quarter)

plt.plot(time[:len(avg_rate)],int_rate_new[0:len(avg_rate)])
plt.plot(time[:len(avg_rate)],long_term_rate)
plt.plot(time[:len(avg_rate)],avg_rate)

plt.legend(['Short Rate','2 year rate','Avg. Rate Paid'])
plt.savefig(figure_dir + 'delayed_response.pdf')
