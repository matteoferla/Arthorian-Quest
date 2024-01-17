from .query import QueryArthor

import contextlib
import pandas as pd
import warnings
from rdkit import RDLogger, Chem
from rdkit.Chem import AllChem, Draw
from IPython.display import display
from typing import Dict, Union


def silence_rdkit():
    """
    Silence RDKit warnings
    """
    logger = RDLogger.logger()
    logger.setLevel(RDLogger.CRITICAL)

def show_experiment(query3d: Chem.Mol, experiment_name: str='', save: bool=False):
    """
    Show the 2D structure of the query molecule, experiment_name and its SMARTS pattern.

    :param query3d:
    :param experiment_name:
    :return:
    """
    query = Chem.Mol(query3d)
    if not experiment_name:
        experiment_name = query.GetProp('experiment') if query.HasProp('experiment') else 'test'
    AllChem.Compute2DCoords(query)
    print(Chem.MolToSmarts(query))
    display(query)
    if save:
        Draw.MolToFile(query, f'{experiment_name}.png')
    return None

