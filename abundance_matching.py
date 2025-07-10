"""
Implementing the abundance matching a la matias.
"""
#TODO: Remove areas. I'm not sure that we need the areas in here. I think the area was just there
# in Matias script because the randoms were for the entire field.

from dataclasses import dataclass
import pandas as pd
import numpy as np
import pylab as plt
from scipy.interpolate import make_splrep

z_bin = np.linspace(0, 0.6, 41)  # 40 bins
z_mid = (z_bin[:-1] + z_bin[1:]) / 2
plot_x = np.linspace(0.001, 0.6, 1000)  # for plotting

LIGHTCONE_PREFIX = "mocks/gama_mock_data/all_lightcones/"

lightcones = [
    "gama0_gals.parquet",
    "gama10_gals.parquet",
    "gama2_gals.parquet",
    "gama4_gals.parquet",
    "gama5_gals.parquet",
    "gama7_gals.parquet",
    "gama3_gals.parquet",
    "gama8_gals.parquet",
    "gama9_gals.parquet",
]


@dataclass
class Region:
    """
    Data structure for gama regions containing area, random file, etc.
    """
    area: float
    random_file_name: str
    ra_range: tuple[float, float]


g09 = Region(area=59.978679332, random_file_name="gama_g09_randoms.txt", ra_range=(129, 141))
g12 = Region(area=59.978679332, random_file_name="gama_g12_randoms.txt", ra_range=(174, 186))
g15 = Region(area=59.978679332, random_file_name="gama_g15_randoms.txt", ra_range=(211.5, 223.5))
g23 = Region(area=50.587431294, random_file_name="gama_g23_randoms.txt", ra_range=(339, 351))

regions = [g09, g12, g15, g23]

GAMA_MAG_LIMIT = 19.65
RANDOM_OVERFACTOR = 400  # The number of times larger the random catalog is.
REDSHIFT_CUT = 0.5


for lightcone in lightcones:
    abundance_matched_fields = []
    for region in regions:
        # Reading in the random data
        df_randoms = pd.read_csv(region.random_file_name)

        # Reading in the gama09 shark galaxies
        df_shark = pd.read_parquet(LIGHTCONE_PREFIX + lightcone)
        df_shark = df_shark[
            (df_shark["ra"] >= region.ra_range[0])
            & (df_shark["ra"] <= region.ra_range[1])
            & (df_shark["total_ap_dust_r_SDSS"] < GAMA_MAG_LIMIT)
            & (df_shark["zobs"] <= REDSHIFT_CUT)
            ]

        df_gama = pd.read_csv(f"gama_galaxy_catalogs/{region.random_file_name.split('_')[1]}_galaxies.dat", sep = ' ')
        # finding the offset in the z-band that is needed to
        nz_random = np.histogram(df_randoms["z"], bins=z_bin, density=False)[0]/ RANDOM_OVERFACTOR
        nz_shark = np.histogram(df_shark["zobs"], bins=z_bin, density=False)[0]
        nz_ratio = nz_random / nz_shark

        mag_difference = np.zeros(len(z_mid))
        for i in range(len(z_bin) - 1):
            z_selection = (df_shark["zobs"] > z_bin[i]) & (df_shark["zobs"] <= z_bin[i + 1])
            number_in_bin = len(np.where(z_selection == True)[0])
            if number_in_bin > 100:  # This seems arbitary
                shark_mags = np.sort(df_shark.loc[z_selection, "total_ap_dust_r_SDSS"].to_numpy(
                        copy=True))
                n_selection = np.sum(shark_mags <= GAMA_MAG_LIMIT)
                idx = int(np.min([np.round(n_selection * nz_ratio[i]), number_in_bin - 1]))
                mag_difference[i] = GAMA_MAG_LIMIT - shark_mags[idx]

        # visualising the offset in each bin.
        mag_diff_spline = make_splrep(z_mid, mag_difference, k=3, s=1)

        plt.scatter(z_mid, mag_difference)
        plt.plot(plot_x, mag_diff_spline(plot_x))
        plt.show()

        # apply the correction to the shark catalog
        df_shark["total_ap_dust_r_SDSS_matched"] = \
            df_shark["total_ap_dust_r_SDSS"] + mag_diff_spline(df_shark["zobs"])
        new_shark = df_shark[df_shark["total_ap_dust_r_SDSS_matched"] <= GAMA_MAG_LIMIT]

        plt.hist(new_shark["zobs"], bins=z_bin, histtype="step", label="Matched SHARK")
        plt.hist(df_shark["zobs"], bins=z_bin, histtype="step", label="Not matched SHARK")
        plt.hist(df_randoms['z']/400, bins=z_bin, histtype='step', label="Randoms")
        plt.hist(df_gama["Z"], bins=z_bin, label="GAMA galaxies")
        plt.legend()
        plt.ylim(0, 5000)
        plt.show()
        print(len(df_shark), len(new_shark), len(df_gama))
        print(len(df_shark)/len(new_shark))
        print(len(df_shark)/len(df_gama))
        print(len(new_shark)/len(df_gama))
        print()
        abundance_matched_fields.append(new_shark)

    abundanced_matched_lightcone = pd.concat(abundance_matched_fields)
    new_name = lightcone.split(".parquet", maxsplit=1)[0] + "_matched.parquet"
    abundanced_matched_lightcone.to_parquet(LIGHTCONE_PREFIX + new_name)
