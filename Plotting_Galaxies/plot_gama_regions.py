"""
Script to simply plot the gama regions to ensure they are correct.
"""

from dataclasses import dataclass
from astropy.io import fits
from astropy.table import Table
import numpy as np
import pylab as plt
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


def polar_to_cartesian(
    radians: float | np.ndarray[float], redshift: float | np.ndarray[float]
) -> tuple[float | np.ndarray[float], float | np.ndarray[float]]:
    """
    Converts the RA and redshift into cartesian coordinates that can be plotted.
    """
    x = redshift * np.cos(radians)
    y = redshift * np.sin(radians)
    return x, y


def read_in_fits_data(fits_file_name: str) -> pd.DataFrame:
    """
    Opens the fits file and reads in the data returning it as a pandas data-frame.
    """
    with fits.open(fits_file_name) as fits_data:
        table = Table(fits_data[1].data)
    df = table.to_pandas()
    return df[df["z"] > 0]  # remove these -9 things


def cut_region(data: pd.DataFrame, region: Region) -> pd.DataFrame:
    """
    Cuts the gama regions out based on the lower_ra an upper_ra. These are set from the
    GAMA website: (https://www.gama-survey.org/dr4/)
    """
    return data[(data["RAcen"] > region.ra_lower) & (data["RAcen"] < region.ra_upper)]


def draw_vertical_line(
    radian_value: float, r_limits: tuple[float, float], **kwargs
) -> None:
    """
    Plots a line polar space but in cartesian coords.
    """
    rs = np.linspace(r_limits[0], r_limits[1], 1000)
    radians = radian_value * np.ones(len(rs))
    line_x, line_y = polar_to_cartesian(radians, rs)
    plt.plot(line_x, line_y, **kwargs)


def draw_horizontal_line(
    radian_range: tuple[float, float], r_value: float, **kwargs
) -> None:
    """
    Draws a horizontal line in polar coordinates.
    """
    plotting_phi = np.linspace(radian_range[0], radian_range[1], 1000)
    line = np.ones(len(plotting_phi)) * r_value
    line_x, line_y = polar_to_cartesian(plotting_phi, line)
    plt.plot(line_x, line_y, **kwargs)



class RegionScatterPlot:
    """
    Main class for plotting wedge diagrams which are scatter points.
    """

    def __init__(
        self,
        ra: np.ndarray,
        redshift: np.ndarray,
        r_max: float = None,
        fig_size: tuple[float, float] = (5, 10),
        **kwargs,
    ) -> None:
        if r_max is not None:
            ra = ra[redshift < r_max]
            redshift = redshift[redshift < r_max]

        phi = np.deg2rad(ra)
        self.alpha = (np.pi / 2) - np.mean(phi)
        self.radians = phi + alpha
        self.r = redshift
        self.r_lim = np.max(self.r)
        self.radian_limits = (np.min(self.radians), np.max(self.radians))

        fig = plt.figure(figsize=fig_size)
        ax = fig.add_subplot()
        ax.set_axis_off()
        x, y = polar_to_cartesian(self.radians, self.r)
        ax.scatter(x, y, **kwargs)

    def plot_border(self, **kwargs) -> None:
        """
        Creates the border around the wedge diagram
        """
        draw_vertical_line(self.radian_limits[0], (0, self.r_lim), **kwargs)
        draw_vertical_line(self.radian_limits[1], (0, self.r_lim), **kwargs)
        draw_horizontal_line(self.radian_limits, self.r_lim, **kwargs)

    def plot_grid(
        self,
        radian_positions: np.ndarray = None,
        r_positions: np.ndarray = None,
        labels: bool = True,
        **kwargs,
    ) -> None:
        """
        Adds a grid. If no values are actually given then the grid is estimated by splitting the
        axis equally into 5 parts.
        """
        if radian_positions is None:
            radian_positions = np.linspace(*self.radian_limits, 7)[1:-1]

        if r_positions is None:
            r_positions = np.linspace(0, self.r_lim, 7)[1:-1]

        for radian_pos in radian_positions:
            draw_vertical_line(radian_pos, (0, self.r_lim), **kwargs)

        for r_pos in r_positions:
            draw_horizontal_line(self.radian_limits, r_pos, **kwargs)

        if labels:
            print(r_positions)
            for rad_pos, r_pos in zip(radian_positions, r_positions):
                text_x, text_y = polar_to_cartesian(self.radian_limits[1], r_pos)
                plt.text(text_x - 0.01, text_y, f"{r_pos:.1f}")

    def labels(self, ra_values: np.ndarray, z_values: np.ndarray, **kwargs) -> None:
        """
        Will place the given labels.
        """
        #TODO: implement the labels
        rad_values = np.deg2rad(ra_values) + self.alpha
        z_axis_points = np.ones(len(rad_values))*self.radian_limits[1]

if __name__ == "__main__":
    INFILE = "../GAMA_galaxies.fits"
    data = read_in_fits_data(INFILE)
    g_09 = cut_region(data, g_09_footprint)
    g_12 = cut_region(data, g_12_footprint)
    g_15 = cut_region(data, g_15_footprint)
    g_23 = cut_region(data, g_23_footprint)

    scatter = RegionScatterPlot(g_09["RAcen"], g_09["z"], 0.5, s=0.1, color="k")
    scatter.plot_border(color="k", lw=3)
    scatter.plot_grid(r_positions=np.arange(0, 0.6, 0.1), color="r", lw=1)
    plt.show()
