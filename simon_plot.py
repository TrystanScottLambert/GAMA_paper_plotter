"""
Script to make Simons plot, comparing the percentage of the total mass of the halo, contained within
the central galaxy as a function of halo mass.
"""
import numpy as np
import seaborn as sns
import pandas as pd
from astropy.table import Table
import pylab as plt

from plotting import start_plot, end_plot

INFILE_GAL = '~/Desktop/my_tools/make_gama_dmu/G3CGal.fits'
INFILE_GROUP = "~/Desktop/my_tools/make_gama_dmu/G3CFoFGroup.fits"


tbl = Table.read(INFILE_GAL)
df_gal = tbl.to_pandas()

tbl = Table.read(INFILE_GROUP)
df_group = tbl.to_pandas()
df_group = df_group[df_group['Nfof'] > 3]

df_gal['lum'] = 10**(-0.4 * df_gal['AbsoluteMagR'])
bcg_ids = df_group['IterCenUberID']

bcg_dfs = df_gal[df_gal['UberID'].isin(bcg_ids)]

df_masses = pd.read_csv('gama_stellar_masses.csv')
df_masses['UberID'] = df_masses['uberID']
bcg_dfs = bcg_dfs.merge(df_masses, on='UberID', how='left')


fraction = np.array(bcg_dfs['StellarMass_50'])/np.array(df_group['MassA'])
plt.scatter(np.log10(df_group['MassAfunc']), np.log10(fraction))
#plt.ylim(0, 1)
plt.xlim(left=11)
plt.show()

plt.hist2d(np.log10(df_group['MassAfunc']), np.log10(fraction), bins=100, cmap='plasma')
plt.xlim(left=11)
plt.show()


start_plot("M_halo", "fraction bcg/halo")
plt.hexbin(np.log10(df_group['MassAfunc']), np.log10(fraction), gridsize=50, cmap='inferno', mincnt=1)
end_plot("simon_plot.png")
plt.close()
