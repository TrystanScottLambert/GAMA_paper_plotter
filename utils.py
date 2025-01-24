"""
Utility functions
"""

from dataclasses import dataclass
import numpy as np
import pandas as pd


@dataclass
class Region:
    """
    GAMA geometric regions (ra/dec limits)
    """

    ra_limits: tuple
    dec_limits: tuple

    def __post_init__(self):
        self.ra_lower = np.min(self.ra_limits)
        self.ra_upper = np.max(self.ra_limits)
        self.dec_lower = np.min(self.dec_limits)
        self.dec_upper = np.max(self.dec_limits)


g_09_footprint = Region(ra_limits=(129.0, 141.0), dec_limits=(-2, 3))
g_12_footprint = Region(ra_limits=(174.0, 186.0), dec_limits=(-3, 2))
g_15_footprint = Region(ra_limits=(211.5, 223.5), dec_limits=(-2, 3))
g_23_footprint = Region(ra_limits=(339.0, 351.0), dec_limits=(-35, 30))

def cut_region(data: pd.DataFrame, region: Region, ra_label: str = 'RAcen') -> pd.DataFrame:
    """
    Cuts the gama regions out based on the lower_ra an upper_ra. These are set from the
    GAMA website: (https://www.gama-survey.org/dr4/)
    """
    return data[(data[ra_label] > region.ra_lower) & (data[ra_label] < region.ra_upper)]
