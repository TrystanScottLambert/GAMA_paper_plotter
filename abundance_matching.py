"""
Implementing the abundance matching a la matias.
"""

import pandas as pd
import numpy as np
import pylab as plt
from scipy.interpolate import make_splrep

GAMA_AREA = 59.978679332 # deg^2 G09 area for testing. TODO: make this for all GAMA regions
GAMA_MAG_LIMIT = 19.65
RANDOM_OVERFACTOR = 400 # The number of times larger the random catalog is.

# reading in the actual gama09 region data
df_gama = pd.read_csv("/Users/00115372/Desktop/GAMA_paper_plotter/gama_galaxy_catalogs/g09_galaxies.dat", sep=" ")

# Reading in the gama09 random data
df_randoms = pd.read_csv("gama_g09_randoms.txt")

# Reading in the gama09 shark galaxies
df_shark  = pd.read_parquet("/Users/00115372/Desktop/GAMA_paper_plotter/mocks/gama_mock_data/all_lightcones/gama0_gals.parquet")
df_shark  = df_shark[(df_shark["ra"]>129) & (df_shark["dec"] < 141)] # selecting the g09 region
df_shark = df_shark[df_shark["total_ap_dust_r_SDSS"] < GAMA_MAG_LIMIT] # applying GAMA cut


# finding the offset in the z-band that is needed to
z_bin = np.linspace(0, 0.6, 41) # 40 bins
z_mid = (z_bin[:-1] + z_bin[1:])/2
plot_x = np.linspace(0.001, 0.6, 1000) # for plotting

nz_random = np.histogram(df_randoms["z"], bins=z_bin, density=False)[0]/RANDOM_OVERFACTOR
nz_shark = np.histogram(df_shark["zobs"], bins=z_bin, density=False)[0]
nz_ratio = nz_random/nz_shark

mag_difference = np.zeros(len(z_mid))
for i in range(len(z_bin) - 1):
    z_selection = (df_shark["zobs"] > z_bin[i]) & (df_shark["zobs"] <= z_bin[i+1])
    number_in_bin = len(np.where(z_selection == True)[0])
    if number_in_bin > 100: # This seems arbitary
        shark_mags = np.sort(df_shark.loc[z_selection, "total_ap_dust_r_SDSS"].to_numpy(copy=True))
        n_selection = np.sum(shark_mags <= GAMA_MAG_LIMIT)
        idx = int(np.min([np.round(n_selection * nz_ratio[i]), number_in_bin -1]))
        mag_difference[i] = GAMA_MAG_LIMIT - shark_mags[idx]

# visualising the offset in each bin.
mag_diff_spline = make_splrep(z_mid, mag_difference, k=3, s=0.1)

plt.scatter(z_mid, mag_difference)
plt.plot(plot_x, mag_diff_spline(plot_x))
plt.show()


# apply the correction to the shark catalog
df_shark["total_ap_dust_r_SDSS_matched"] = df_shark["total_ap_dust_r_SDSS"] + mag_diff_spline(df_shark["zobs"])
new_shark = df_shark[df_shark["total_ap_dust_r_SDSS_matched"] <= GAMA_MAG_LIMIT]

plt.hist(new_shark["zobs"], bins=z_bin, histtype='step', label='Matched SHARK')
plt.hist(df_shark["zobs"], bins=z_bin, histtype='step', label= 'Not matched SHARK')
#plt.hist(df_randoms['z']/400, bins=z_bin, histtype='step', label="Randoms")
plt.hist(df_gama["Z"], bins=z_bin, label="GAMA galaxies")
plt.legend()
plt.ylim(0, 5000)
plt.show()

