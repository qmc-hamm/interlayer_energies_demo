from .generate_geometry import create_graphene_geom
import pkg_resources
import pandas as pd


def load_data():
    """Return a dataframe 

    Contains the following fields:
        energy: eV/atom
        energy_err: eV/atom
        finite_size_correction : eV/atom : the difference between the largest calculation done and the extrapolated infinite value
        d: distance between the layers (angstrom)
        disregistry:
            0.0 : AB stacking
            1/6: SP stacking
            2/3: AA stacking

    """
    # This is a stream-like object. If you want the actual info, call
    # stream.read()
    stream = pkg_resources.resource_stream(__name__, 'data/qmc.csv')
    return pd.read_csv(stream)
