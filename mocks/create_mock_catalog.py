"""
Script that converts the parquet files that we have downloaded from SHARK on pawsey into a format
that is more accepted by the FoFR group finder. In particular renaming the groupids to be 1..N and
renaming the column names to the keywords that are looked for in the R code.
"""

import pandas as pd
import numpy as np

# from astropy.table import Table
def assign_g_field(ra: float) -> str:
    """
    Assigns which gama field the ra value is in.
    """
    if 128 < ra < 142:
        field = 'g09'
    elif 173 < ra < 187:
        field = 'g12'
    elif 211 < ra < 224:
        field = 'g15'
    elif 338 < ra < 352:
        field = 'g23'
    else:
        field = 'NOT IN GAMA'
    return field

def rename_groups(
    data_frame: pd.DataFrame, group_id_col_name: str, static_id: int
) -> pd.DataFrame:
    """
    Renaming the data_frame groups from 1 .. N instead of the stingray naming convention.
    This is the way that the FoFR finder wants to work.
    """
    unique_ids = data_frame[group_id_col_name].unique()
    isolated_mask = data_frame[group_id_col_name] == static_id
    grouped_ids = sorted(unique_ids[unique_ids != static_id])

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

def rename_groups_new(data_frame: pd.DataFrame, group_id_col_name: str) -> pd.DataFrame:
    """
    Renames group IDs in the DataFrame:
    - Keeps isolated galaxies (group_id == -1) unchanged.
    - Renumbers all other group IDs to be sequential from 1 to N.
    """
    # Copy to avoid modifying original
    df = data_frame.copy()

    # Identify unique group IDs excluding -1
    valid_ids = sorted(df.loc[df[group_id_col_name] != -1, group_id_col_name].unique())

    # Create new mapping: old_id -> new_id (starting from 1)
    new_id_mapping = {old_id: new_id for new_id, old_id in enumerate(valid_ids, start=1)}

    # Apply mapping, leave -1 untouched
    df["id_group_sky"] = df[group_id_col_name].map(new_id_mapping).fillna(-1).astype(int)

    return df


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
    outfile = "gama_gals_for_R.parquet"
    lightcones_number = [0, 2, 3, 4, 5, 7, 8, 9, 10]

    lightcones = []
    for lightcone_number in lightcones_number:
        print(f"Now doing lightcone: {lightcone_number}")
        infile = f"gama_mock_data/all_lightcones/gama{lightcone_number}_gals.parquet"
        df = pd.read_parquet(infile)
        df = rename_groups_new(df, "id_group_sky")
        df = rename_ids_col_names(df)
        df = df[df["total_ap_dust_r_SDSS_matched"] < 19.65]
        df = df[df["zobs"] < 0.5]
        df = df[df["ra"] > 128]  # removing the gama02 region
        df["LC"] = np.ones(len(df)).astype(int) * lightcone_number
        lightcones.append(df)

    df_all = pd.concat(lightcones)
    df_all['g_field'] = df_all['ra'].apply(assign_g_field)
    df_all['lightcone_gamafield'] = df_all['LC'].astype(str) + df_all['g_field']
    df_all.to_parquet(outfile, index=False)

if __name__ == "__main__":
    main()
