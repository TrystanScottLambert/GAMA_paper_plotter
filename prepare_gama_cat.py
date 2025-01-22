"""
Preparing the GAMA catalog from its fits file to it's .dat file. This is 

This is a *script* it's only meant to put the GAMA data into the .dat file so the current 
implementation of the R code will read it. 
"""

from astropy.io import fits
from astropy.table import Table

from utils import cut_region, g_09_footprint

def main():
    """
    main script
    """
    infile = 'GAMA_galaxies.fits'
    with fits.open(infile) as fits_data:
        data = Table(fits_data[1].data)
        df = data.to_pandas()

    df = df[(df['z'] > 0.01) & (df['z'] < 1)] # removing weird negative velocities.
    df = cut_region(df, g_09_footprint)  # we are only looking at g09 for now.
    df = df[['RAcen', 'Deccen', 'z', 'mag']]
    df = df[df['mag'] < 19.65]
    df.rename(columns={'RAcen':"RA", 'Deccen': 'DEC', 'z':'Z', 'mag':'Rpetro'}, inplace=True)
    df.to_csv('GAMA_galaxies.dat', sep=' ', index=False) #full gama catalog

    # creating a test case for g09. It seems that the we have to go region by region for the 
    # R group finder

if __name__ == '__main__':
    main()
