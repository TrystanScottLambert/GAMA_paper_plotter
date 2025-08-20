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
df_group = df_group[df_group['Nfof'] > 5]
df_group = df_group[df_group['Zfof'] < 0.2]

bcg_ids = df_group['BCGUberID']

bcg_dfs = df_gal[df_gal['UberID'].isin(bcg_ids)]

df_masses = pd.read_csv('gama_stellar_masses.csv')
df_masses['UberID'] = df_masses['uberID']
bcg_dfs = bcg_dfs.merge(df_masses, on='UberID', how='left')


fraction = np.array(bcg_dfs['StellarMass_50'])/np.array(df_group['MassA'])
start_plot("M_halo", "fraction bcg/halo")
plt.hexbin(np.log10(df_group['MassA']), fraction, gridsize=50, cmap='inferno', mincnt=1)
#plt.xlim(8, 16)
#plt.ylim(-6, 2)
end_plot("simon_plot.png")
plt.close()


"""
#Trying the same analysis except using the old data.
"""

tbl = Table.read("/Users/00115372/Desktop/GAMA_paper_plotter/old_gama_data/ProSpectv03.fits")
df_prospect = tbl.to_pandas()

tbl = Table.read("/Users/00115372/Desktop/GAMA_paper_plotter/old_gama_data/G3CFoFGroupv10.fits")
df_groups_10 = tbl.to_pandas()
df_groups_10 = df_groups_10[df_groups_10['Nfof'] > 5]
df_groups_10 = df_groups_10[df_groups_10['Zfof'] < 0.2]

bcg_ids = np.array(df_groups_10['BCGCATAID'])

all_ids = np.array(df_prospect['CATAID'])
idx = []
bad_group_idx = []
for i, _id in enumerate(bcg_ids):
    _idx = np.where(all_ids == _id)[0]
    if len(_idx) != 0:
        idx.append(_idx[0])
    else:
        bad_group_idx.append(i)
bad_group_idx = np.array(bad_group_idx)
good_group_idx = np.setdiff1d(np.arange(len(df_groups_10)), bad_group_idx)
idx = np.array(idx)

bad_groups = df_groups_10.iloc[bad_group_idx]
df_groups_10 = df_groups_10.iloc[good_group_idx]

stellar_mass_bcg = np.array(df_prospect.iloc[idx]['StellarMass_50'])
halo_mass_a = np.array(df_groups_10['MassA'])
halo_mass_afunc = np.array(df_groups_10['MassAfunc'])

start_plot("Halo Mass [Msol]", "Fraction of masss in BCG [percent]")
plt.hexbin(np.log10(halo_mass_a), stellar_mass_bcg/halo_mass_a, gridsize=50, cmap='inferno', mincnt=1)
#plt.xlim(8, 16)
#plt.ylim(-6, 2)
end_plot("simon_plot_old_group_cat.png")
plt.close()


"""
Working out the same plot but for Shark
"""

#INFILE_SHARK = "/Users/00115372/Desktop/mock_catalogs/offical_waves_mocks/v0.3.1/waves_wide_gals.parquet"
#df_shark = pd.read_parquet(INFILE_SHARK)


# TODO: Shark needs the Type included for me to work this trend out directly from simulations.
# TODO: Once shark has worked that out then we can run Nesie on Shark and see how that differs.
