"""
Interrogating the group data.
"""

import pandas as pd
import pylab as plt
import numpy as np
import pyvista as pv
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM

from Plotting_Galaxies.plot_gama_regions import RegionScatterPlot
from utils import cut_region, g_23_footprint


def quick_distribution(data_frame: pd.DataFrame, col_name: str, bins: int = 100, **kwargs) -> None:
    """
    Quickly plots the distribution for checking.
    """
    plt.hist(data_frame[col_name], bins=bins, **kwargs)
    plt.show()

def add_xyz(data_frame: pd.DataFrame, ra_col: str, dec_col: str, z_col: str) -> pd.DataFrame:
    """
    Calculates the x, y, z coordinates for the data frame.
    """
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    distance = cosmo.comoving_distance(np.array(data_frame[z_col]))
    c = SkyCoord(
        ra=np.array(data_frame[ra_col])*u.degree,
        dec=np.array(data_frame[dec_col])*u.degree,
        distance=distance
        )
    data_frame['x'] = c.cartesian.x.value
    data_frame['y'] = c.cartesian.y.value
    data_frame['z'] = c.cartesian.z.value
    return data_frame


if __name__ == '__main__':
    gal_ids, group_ids = np.loadtxt('g09_galaxy_linking_table.csv', delimiter=',', unpack=True, skiprows=1)

    df_galaxies = pd.read_csv('gama_galaxy_catalogs/g09_galaxies.dat', sep='\s+')
    add_xyz(df_galaxies, 'RA', 'DEC', 'Z')
    df_group_galaxies = df_galaxies.iloc[gal_ids - 1]


    infile = 'g23_group_catalog.csv'
    df = pd.read_csv(infile)
    add_xyz(df, 'IterCenRA', 'IterCenDEC', 'MedianZ')


    mock_gal_ids, mock_group_ids = np.loadtxt('testing_broken_tuning_linking_table.csv')
    df_mock_galaxies = pd.read_parquet('mocks/gama_gals_for_R.parquet')
    add_xyz(df_mock_galaxies, 'ra', 'dec', 'zobs')
    df_mock_group_galaxies = df_mock_galaxies[mock_gal_ids -1]

    

    scatter = RegionScatterPlot(df['IterCenRA'], df['MedianZ'], 1, s=np.log10(df['Mult'])*10, alpha=0.5, facecolor='none', edgecolors='k')
    scatter.plot_border(color='k', lw=3)
    scatter.plot_grid(color='r', alpha=0.1)
    plt.savefig('current_groups.png')

    # Read in the GAMA thing. 
    df_aaron = pd.read_csv('GAMA_groups_aaron.csv')
    df_aaron_g09 = cut_region(df_aaron, g_23_footprint, ra_label='IterCenRA')
    scatter = RegionScatterPlot(df_aaron_g09['IterCenRA'], df_aaron_g09['IterCenZ'], s=np.log10(df_aaron_g09['Nfof'])*10, alpha=0.5)
    scatter.plot_border(color='k', lw=3)
    scatter.plot_grid(color='r', alpha=0.1)
    plt.savefig('aaron_groups.png')

    THREE_D_PLOT = False
    if THREE_D_PLOT:
        plotter = pv.Plotter()

        # Loop through each row in the dataframe and add a sphere
        for _, row in df.iterrows():
            center = (row['x'], row['y'], row['z'])  # Center of the sphere
            radius = row['Rad50']  # Radius of the sphere
            sphere = pv.Sphere(radius=radius, center=center)  # Create a sphere
            plotter.add_mesh(sphere, style='wireframe', color='r', opacity=0.2)

        galaxy_points = df_group_galaxies[['x', 'y', 'z']].to_numpy()  # Convert to NumPy array
        point_cloud = pv.PolyData(galaxy_points)
        plotter.add_mesh(point_cloud, color='black', point_size=2, render_points_as_spheres=True)

        # Add labels or enhance visualization (optional)
        #plotter.add_axes()  # Add coordinate axes
        #plotter.show_bounds(grid='back')  # Show bounding box and grid

        # Render the plot
        plotter.show()
