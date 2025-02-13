"""
Testing the create_mock_catalog.py
"""

from io import BytesIO
import unittest
import pandas as pd

from create_mock_catalog import rename_ids_col_names, rename_groups


class TestRenameFunctions(unittest.TestCase):
    """
    Main testing unit for the create_mock_catalog script.
    """

    def setUp(self):
        # Create a test dataframe
        self.df = pd.DataFrame(
            {
                "id_group_sky": [
                    -1,
                    -1,
                    -1,
                    6292300000001,
                    6292300000001,
                    6292300000001,
                ],
                "id_galaxy_sky": [101, 102, 103, 104, 105, 106],
            }
        )

        # Save to a Parquet file in memory
        self.buffer = BytesIO()
        self.df.to_parquet(self.buffer)
        self.buffer.seek(0)

    def test_rename_group_ids(self):
        """
        Testing that we rename the group ids correctly.
        """
        df_test = pd.read_parquet(self.buffer)
        df_transformed = rename_groups(df_test)
        expected_group_ids = [1, 2, 3, 4]
        self.assertListEqual(df_transformed["GroupID"].unique().tolist(), expected_group_ids)

    def test_rename_column_names(self):
        """
        Testing that the colums are renamed correctly.
        """
        df_test = pd.read_parquet(self.buffer)
        df_transformed = rename_ids_col_names(df_test)
        new_names = df_transformed.columns.to_list()
        correct = ["GroupID", "CATAID"]
        for cor, new in zip(correct, new_names):
            self.assertEqual(cor, new)


if __name__ == "__main__":
    unittest.main()
