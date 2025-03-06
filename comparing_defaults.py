"""
Looking at what the default parameters look like on SHARK mock versus the current 
best fit parameters.
"""
import pandas as pd
import pylab as plt
import numpy as np

best_fit_group_catalog = pd.read_csv('best_testing_group_catalog.csv')
best_fit_linking_table = pd.read_csv('best_testing_galaxy_linking_table.csv')

default_fit_group_catalog = pd.read_csv('default_testing_group_catalog.csv')
default_fit_linking_table = pd.read_csv('default_testing_galaxy_linking_table.csv')


plt.hist(best_fit_group_catalog['MedianZ'], histtype='step', bins=np.arange(0, 0.5, 0.01), label='Best Fit')
plt.hist(default_fit_group_catalog['MedianZ'], histtype='step', bins=np.arange(0, 0.5, 0.01), label='Default')
plt.legend()
plt.show()
plt.hist(np.log10(best_fit_group_catalog['Mass']), histtype='step', bins=np.arange(6, 15, 0.2), label='Best Fit')
plt.hist(np.log10(default_fit_group_catalog['Mass']), histtype='step', bins=np.arange(6, 15, 0.2), label='Default')
plt.legend()
plt.show()
plt.hist(best_fit_group_catalog['Mult'], histtype='step', bins=np.arange(20), label='Best Fit')
plt.hist(default_fit_group_catalog['Mult'], histtype='step', bins=np.arange(20), label='Default')
plt.legend()
plt.yscale('log')
plt.show()
