"""
Combines the random catalogs that were generated for the four different fields
"""

import pandas as pd

infiles = [
    "gama_g09_randoms.txt",
    "gama_g12_randoms.txt",
    "gama_g15_randoms.txt",
    "gama_g23_randoms.txt",
    ]

dfs = [pd.read_csv(file) for file in infiles]
combined_df = pd.concat(dfs)
combined_df.to_csv('gama_combined_randoms.txt', sep=' ', index=False)
