"""
Inpsecting that the mock group catalog looks like a catalog.
"""

import pandas as pd
import numpy as np
import pylab as plt
from astropy.table import Table


def prepare_shark_lightcone(infile: str, outfile: str) -> None:
    """
    Reads in the given shark lightcone, applies survey limits, renames groups so that all isolated
    galaxies have -1 and the group-ids are in sequential order. Ignoring the i64 issue.
    """
    df = pd.read_parquet(infile)
    df = df[df['ra'] > 100]
    df = df[df['zobs'] < 0.5]
    df = df[df['total_ap_dust_r_SDSS'] < 19.8]

    unique_group_ids, counts = np.unique(df['id_group_sky'], return_counts=True)
    unique_group_ids, counts = unique_group_ids[1:], counts[1:]

    isolated_ids = unique_group_ids[np.where(counts == 1)]
    group_ids = unique_group_ids[np.where(counts>1)]

    # Preparing shark for tuning
    df.loc[df['id_group_sky'].isin(isolated_ids), 'id_group_sky'] = -1
    new_group_ids = np.arange(len(group_ids)) + 1
    id_mapping = dict(zip(group_ids, new_group_ids))
    df['id_group_sky'] = df['id_group_sky'].replace(id_mapping)

    #df = df[df['id_group_sky'] != -1]

    df.to_parquet(outfile)

def prepare_shark_lightcone_new(infile: str, outfile: str) -> None:
    """
    Trying to fix bug somewhere with tuning being much lower for shark.
    """
    df = pd.read_parquet(infile)
    df = df[df['ra'] > 100]
    df = df[df['zobs'] < 0.5]
    df = df[df['total_ap_dust_r_SDSS'] < 19.8]
    print(len(df))
    #df = df.sample(n=100_000)
    df = df[df['id_group_sky'] != -1]
    print(len(df))

    grouped_ids = sorted(df.loc[df['id_group_sky'] != -1, 'id_group_sky'].unique())

    new_id_mapping = {old_id: new_id for new_id, old_id in enumerate(grouped_ids, start=1)}

    df['id_group_sky'] = df['id_group_sky'].map(lambda x: new_id_mapping.get(x, x))
    df.to_parquet(outfile)

def rename_groups(data_frame: pd.DataFrame) -> pd.DataFrame:
    """
    Renaming the data_frame groups from 1 .. N instead of the stingray naming convention.
    This is the way that the FoFR finder wants to work.
    """
    unique_ids = data_frame["id_group_sky"].unique()
    isolated_mask = data_frame["id_group_sky"] == -1
    grouped_ids = sorted(unique_ids[unique_ids != -1])

    number_isolated_galaxies = isolated_mask.sum()
    data_frame.loc[isolated_mask, "id_group_sky"] = range(
        1, number_isolated_galaxies + 1
    )

    new_id_mapping = {
        old_id: new_id
        for new_id, old_id in enumerate(grouped_ids, start=number_isolated_galaxies + 1)
    }
    data_frame["id_group_sky"] = data_frame["id_group_sky"].replace(new_id_mapping)
    return data_frame

def prepare_shark_dumb(infile: str, outfile: str) -> None:
    """
    i ahte my life
    """
    df = pd.read_parquet(infile)
    df = rename_groups(df)
    df = df[df['ra'] > 100]
    df = df[df['zobs'] < 0.5]
    df = df[df['total_ap_dust_r_SDSS'] < 19.8]
    df.to_parquet(outfile)



def prepare_galform_lightcone(infile: str, outfile: str) -> None:
    """
    Renames the single groups to -1 and writes to parquet.
    """
    galform_data = Table.read(infile)
    galform_data = galform_data.to_pandas()
    galform_data['GroupID'][np.where(galform_data['GroupID'] == 0)] = -1
    galform_data.to_parquet(outfile)


if __name__ == '__main__':
    #Prepare Shark
    lightcones = [0, 2, 3, 4, 5, 7, 8, 9, 10]
    infiles = [f'gama_mock_data/all_lightcones/gama{num}_gals.parquet' for num in lightcones]
    outfiles = [f'shark_for_R_{num}.parquet' for num in lightcones]
    for _in, out in zip(infiles, outfiles):
        print('Doing: ', _in)
        prepare_shark_lightcone(_in, out)

    ##Prepare GALFORM
    #prepare_galform_lightcone('G3CMockGALv04.fits', 'galform_gals_for_R.parquet')

    '''df_shark = pd.read_parquet('shark_for_R_0.parquet')
    df_galform = pd.read_parquet('galform_gals_for_R.parquet')
    df_shark_0g09 = df_shark[df_shark['ra'] < 150]
    df_galform_0g09 = df_galform[(df_galform['RA'] < 150) & (df_galform['Volume'] == 1)]

    # n(z) distribution
    z_bins = np.arange(0, 0.6, 0.01)
    plt.hist(df_shark_0g09['zobs'], bins=z_bins, histtype='step', label='Shark')
    plt.hist(df_galform_0g09['Zspec'], bins=z_bins, histtype='step', lw=3, label='Galform')
    plt.legend()
    plt.show()

    # n(z) centrals and satelites in shark. 
    centrals = df_shark[df_shark['type'] == 0]
    satelites = df_shark[df_shark['type'] == 1]
    plt.hist(centrals['zobs'], bins=z_bins, histtype='step', label='centrals')
    plt.hist(satelites['zobs'], bins=z_bins, histtype='step', lw=3, label='satelites')
    plt.legend()
    plt.show()


    # multiplicity distribution
    shark_group_ids, shark_group_counts = np.unique(df_shark_0g09['id_group_sky'], return_counts=True)
    shark_group_ids, shark_group_counts =  shark_group_ids[1:], shark_group_counts[1:]

    galform_group_ids, galform_group_counts = np.unique(df_galform_0g09['GroupID'], return_counts=True)
    galform_group_ids, galform_group_counts = galform_group_ids[1:], galform_group_counts[1:]

    print('The number of unique groups in Shark: ', len(np.where(shark_group_counts > 5)[0]))
    print('The number of unique groups in Galform: ', len(np.where(galform_group_counts > 5)[0]))

    count_bins = np.arange(2, 50, 1)
    plt.hist(shark_group_counts, bins=count_bins, histtype='step', label='Shark')
    plt.hist(galform_group_counts, bins=count_bins, histtype='step', label='Galform')
    plt.yscale('log')
    plt.legend()
    plt.show()

    # side by side comparison of ra and dec

    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.scatter(df_galform_0g09['RA'], df_galform_0g09['DEC'], s=0.1)

    ax = fig.add_subplot(212)
    ax.scatter(df_shark_0g09['ra'], df_shark_0g09['dec'], s=0.1)
    plt.show()

    shark_group_ids = shark_group_ids[np.where(shark_group_counts>2)]
    galform_group_ids = galform_group_ids[np.where(galform_group_counts>2)]
    shark_groups = df_shark_0g09[df_shark_0g09['id_group_sky'].isin(shark_group_ids)]
    galform_groups = df_galform_0g09[df_galform_0g09['GroupID'].isin(galform_group_ids)]

    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.scatter(galform_groups['RA'], galform_groups['DEC'], s=0.1)

    ax = fig.add_subplot(212)
    ax.scatter(shark_groups['ra'], shark_groups['dec'], s=0.1)
    plt.show()'''
