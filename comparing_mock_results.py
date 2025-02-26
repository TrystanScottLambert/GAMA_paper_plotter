"""
Comparing the results of running on the mock catalog and investingating that the mock catalog is
actually giving us a reasonable mock catalog.
"""

import pandas as pd
import pylab as plt
import numpy as np


if __name__ == '__main__':
    # first compare that the mock actually looks like the gama field.
    mock_df = pd.read_parquet('mocks/gama_gals_for_R.parquet')
    actual_data = pd.read_csv('gama_galaxy_catalogs/g09_galaxies.dat', sep=' ')
    actual_linking_table = pd.read_csv('g09_galaxy_linking_table.csv')
    actual_groups = pd.read_csv('g09_group_catalog.csv')

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


