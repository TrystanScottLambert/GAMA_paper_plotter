"""
Preparing the GAMA catalog from its fits file to it's .dat file. This is 

This is a *script* it's only meant to put the GAMA data into the .dat file so the current 
implementation of the R code will read it. 
"""

from astropy.io import fits
from astropy.table import Table

def main():
    """
    main script
    """
    infile = 'GAMA_galaxies.fits'
    with fits.open(infile) as fits_data:
        data = Table(fits_data[1].data)
        df = data.to_pandas()
    
    df = df[(df['z'] > 0.01) & (df['z'] < 1)] # removing weird negative velocities.
    df = df[['RAcen', 'Deccen', 'z', 'mag']]
    df = df.rename(columns={'RAcen':"RA", 'Deccen': 'DEC', 'z':'Z', 'mag':'Rpetro'})
    df.to_csv('GAMA_galaxies.dat', sep=' ', index=False)
if __name__ == '__main__':
    main()
