"""
Script that converts the parquet files that we have downloaded from SHARK on pawsey into a format
that is more accepted by the FoFR group finder. In particular renaming the groupids to be 1..N and
renaming the column names to the keywords that are looked for in the R code.
"""

import pandas as pd


def rename_groups(data_frame: pd.DataFrame) -> pd.DataFrame:
    """
    Renaming the data_frame groups from 1 .. N instead of the stingray naming convention.
    This is the way that the FoFR finder wants to work.
    """
    unique_ids = data_frame["id_group_sky"].unique()
    isolated_mask = data_frame["id_group_sky"] == -1  # isolated gals in stingray are -1
    grouped_ids = sorted(unique_ids[unique_ids != -1])

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
    data_frame.rename(columns={"id_group_sky": "GroupID", "id_galaxy_sky": "CATAID"})
    return data_frame


if __name__ == "__main__":
    INFILE = "gama_mock_data/gama_gals.parquet"
    OUTFILE = "gama_gals_for_R.parquet"
    df = pd.read_parquet(INFILE)
    df = rename_groups(df)
    df = rename_ids_col_names(df)
    df.to_parquet(OUTFILE, index=False)
