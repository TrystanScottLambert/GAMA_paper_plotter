"""
Script that converts the parquet files that we have downloaded from SHARK on pawsey into a format
that is more accepted by the FoFR group finder. In particular renaming the groupids to be 1..N and
renaming the column names to the keywords that are looked for in the R code.
"""

import pandas as pd
import pylab as plt
#from astropy.table import Table


def rename_groups(data_frame: pd.DataFrame, group_id_col_name: str, static_id: int) -> pd.DataFrame:
    """
    Renaming the data_frame groups from 1 .. N instead of the stingray naming convention.
    This is the way that the FoFR finder wants to work.
    """
    unique_ids = data_frame[group_id_col_name].unique()
    isolated_mask = data_frame[group_id_col_name] == static_id
    grouped_ids = sorted(unique_ids[unique_ids != static_id])

    number_isolated_galaxies = isolated_mask.sum()
    data_frame.loc[isolated_mask, "id_group_sky"] = range(1, number_isolated_galaxies + 1)

    new_id_mapping = {
        old_id: new_id
        for new_id, old_id in enumerate(grouped_ids, start=number_isolated_galaxies + 1)
    }
    data_frame["id_group_sky"] = data_frame["id_group_sky"].replace(new_id_mapping)
    return data_frame


def rename_ids_col_names(data_frame: pd.DataFrame) -> pd.DataFrame:
    """
    Renames the data_frame columns for galaxy ids and group ids to keep consitent with the FoFR
    naming conventions which FoFempint looks for when comparing to mock catalogs.
    """
    df_new = data_frame.rename(columns={"id_group_sky": "GroupID", "id_galaxy_sky": "CATAID"})
    return df_new

def main():
    """
    Rename the ids and columns
    """
    infile = "gama_mock_data/all_lightcones/gama0_gals_matched.parquet"
    outfile = "gama_gals_for_R.parquet"
    df = pd.read_parquet(infile)
    df = rename_groups(df, 'id_group_sky', -1)
    df = rename_ids_col_names(df)
    df = df[df['total_ap_dust_r_SDSS_matched'] < 19.65]
    df.to_parquet(outfile, index=False)

    # Doing galform 
    #galform_infile = 'G3CMockGalv04.fits'
    #galform_outfile = 'gama_gals_for_R_galform.parquet'
    #df_galform = Table.read(galform_infile).to_pandas()
    #df_galform = df_galform[(df_galform['RA'] < 141) & (df_galform['Rpetro'] < 19.65) & (df_galform['Volume'] == 1)]
    #df_galform_new = df_galform.rename(columns={"GalID": "CATAID"})
    #df_galform_new.to_parquet(galform_outfile)

if __name__ == "__main__":
    main()
