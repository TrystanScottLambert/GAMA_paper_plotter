"""
Comparing the results of running on the mock catalog and investingating that the mock catalog is
actually giving us a reasonable mock catalog.
"""

import pandas as pd
import pylab as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table


if __name__ == '__main__':
    # first compare that the mock actually looks like the gama field.
    mock_df = pd.read_parquet('mocks/gama_gals_for_R.parquet')
    actual_data = pd.read_csv('gama_galaxy_catalogs/g09_galaxies.dat', sep=' ')
    actual_linking_table = pd.read_csv('best_testing_galaxy_linking_table.csv')
    actual_groups = pd.read_csv('best_testing_group_catalog.csv')

    number_of_real_isolated = len(actual_data) - len(actual_linking_table)
    number_of_mock_isolated = len(np.where(np.unique(mock_df['GroupID'], return_counts=True)[1] == 1)[0])
    mock_group_dist = np.unique(mock_df['GroupID'], return_counts=True)[1]
    real_group_dist = np.array(actual_groups['Mult'])
    print(f'Number of real singles: {number_of_real_isolated/len(actual_data)}')
    print(f'Number of mock singles: {number_of_mock_isolated/len(mock_df)}')

    # Check one. The distribution of the multiplicity
    bins = np.arange(3, 50, 1)
    plt.hist(mock_group_dist, bins = bins, histtype='step', lw=3, label='mock')
    plt.hist(real_group_dist, bins=bins, histtype='step', lw=2, label='real')
    plt.yscale('log')
    plt.legend()
    plt.xlabel('Multiplicity', fontsize=20)
    plt.ylabel('Counts', fontsize=20)
    plt.show()

    # Check
    # Check two. The host halo mass for all the isolated galaxies from the parquet files.
    # So we are going to use the exact values of the actual isolated galaxies.
    mock_raw_df = pd.read_parquet('mocks/gama_mock_data/gama_gals.parquet')
    mock_raw_gal_dist = np.unique(mock_raw_df['id_group_sky'], return_counts=1)[1] # including the true isolated.
    true_isolated_df = mock_raw_df[mock_raw_df['id_group_sky'] == -1]

    selected_isolated = len(np.where(mock_raw_gal_dist == 1)[0])
    number_apparent_isolated = selected_isolated + len(true_isolated_df)

    number_weird_isolated = len(np.where(true_isolated_df['mvir_hosthalo'] > 10e12)[0])
    print('The percentage of galaxies identified as isolated is: ', number_apparent_isolated/len(mock_raw_df))
    print('The total number of truly isolated galaxies: ', len(true_isolated_df))
    print('The total number of selected isolated galaxies: ', selected_isolated)
    print('Number of true isolated galaxies greater than 10^12: ', number_weird_isolated)
    plt.hist(np.log10(true_isolated_df['mvir_hosthalo']), bins = np.arange(10, 16, 0.1))
    plt.xlabel('log(host virial Mass [Msol])', fontsize=20)
    plt.ylabel('Counts', fontsize=20)
    plt.show()


    #checking the halo mass distribution. 
    mass_bins = np.arange(10, 15, 0.1)
    plt.hist(np.log10(mock_raw_df[mock_raw_df['mvir_hosthalo'] != -1]['mvir_hosthalo']), bins = mass_bins, histtype='step', label='Not Isolated', density=True)
    plt.hist(np.log10(true_isolated_df['mvir_hosthalo']), bins = mass_bins, histtype='step', label='True isolated', density=True)
    plt.legend()
    plt.xlabel('log10(mvir_host_halo)', fontsize=20)
    plt.ylabel('Counts', fontsize=20)
    plt.show()

    # Comparing the Galform mock catalog to the current mock catalog and real data.
    # Reading in Aarons group results
    df_aaron = pd.read_csv('GAMA_groups_aaron.csv')
    df_aaron = df_aaron[(df_aaron['IterCenRA'] > 100) & (df_aaron['IterCenRA'] < 150)]
    df_default_testing = pd.read_csv('default_testing_group_catalog.csv')
    # Reading in galform
    df_galform = Table.read('mocks/G3CMockGalv04.fits').to_pandas()
    df_galform = df_galform[(df_galform['RA'] < 141) & (df_galform['Rpetro'] < 19.65) & (df_galform['Volume'] == 1)]

    gal_form_counts = np.unique(df_galform['GroupID'], return_counts=True)[1]
    print('Number of ungrouped galaxies are: ', gal_form_counts[0]/len(df_galform))
    bins = np.arange(3, 50, 1)
    plt.hist(mock_group_dist, bins = bins, histtype='step', lw=3, label='SHARK')
    #plt.hist(real_group_dist, bins=bins, histtype='step', lw=2, label='real')
    plt.hist(df_aaron['Nfof'], bins=bins, histtype='step', label='Robotham+2011')
    plt.hist(gal_form_counts[1:], bins=bins, histtype='step', lw=4, label='GalForm')
    plt.hist(df_default_testing['Mult'], bins=bins, histtype='step', label='default settings on mock')
    plt.legend()
    plt.yscale('log')
    plt.xlabel('Multiplicity', fontsize=20)
    plt.ylabel('Counts', fontsize=20)
    plt.show()

    # looking at the n(z) of published, galform, shark, and defaults on shark.
    z_bins = np.arange(0, 0.5, 0.01)
    plt.hist(df_aaron['IterCenZ'], bins =z_bins, histtype='step', label='published')
    plt.hist(df_default_testing['MedianZ'], bins=z_bins, histtype='step', label='default on shark')
    plt.legend()
    plt.show()


    # Testing making the ids ourselves.
    mock_raw_df['our_ids'] = \
        np.array(mock_raw_df['tile']).astype(str) + '-' +\
        np.array(mock_raw_df['subvolume']).astype(str) + '-' + \
        np.array(mock_raw_df['id_halo_sam']).astype(str) + '-' +\
        np.array(mock_raw_df['snapshot']).astype(str)
