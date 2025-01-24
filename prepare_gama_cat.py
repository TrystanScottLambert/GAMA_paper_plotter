"""
Preparing the GAMA catalog from its fits file to it's .dat file. This is 

This is a *script* it's only meant to put the GAMA data into the .dat file so the current 
implementation of the R code will read it. 
"""

from astropy.io import fits
from astropy.table import Table
import numpy as np

from utils import cut_region, g_09_footprint

def convert_jy_to_abmag(flux_in_jansky: np.ndarray) -> np.ndarray:
    """
    Converts the jansky flux from gama to an AB mag 
    see (https://pdn4kd.github.io/2020/10/28/janskyabmag.html)
    """
    return -2.5*np.log10(flux_in_jansky) + 8.90

def main():
    """
    main script
    """
    # SELECT UberID,RAcen,Deccen,flux_rt, Z FROM gkvScienceCatv02
    # WHERE NQ>2 AND SC>6 AND duplicate=0 AND mask=0 AND starmask=0
    infile = 'GAMA_galaxies.fits'
    with fits.open(infile) as fits_data:
        data = Table(fits_data[1].data)
        df = data.to_pandas()

    df = df[(df['Z'] > 0.01) & (df['Z'] < 1)] # velocity range.
    df['mag_r'] = convert_jy_to_abmag(df['flux_rt'])

    df = cut_region(df, g_09_footprint)  # we are only looking at g09 for now.
    df = df[['UberID', 'RAcen', 'Deccen', 'Z', 'mag_r']]
    df.rename(columns={'RAcen':"RA", 'Deccen': 'DEC', 'z':'Z', 'mag_r':'Rpetro'}, inplace=True)
    df.to_csv('GAMA_galaxies.dat', sep=' ', index=False) #full gama catalog

    # creating a test case for g09. It seems that the we have to go region by region for the 
    # R group finder

if __name__ == '__main__':
    main()
