__all__ = ['QueryArthor']

import requests
import warnings
import pandas as pd
from typing import List
from rdkit import Chem
from rdkit.Chem import PandasTools, Draw, AllChem


class QueryArthor:
    """
    Query class for Arthorian Quest.
    It queries arthor.docking.org via ``.retrieve``

    See https://arthor.docking.org/api.html for the API endpoints used

    .. code-block:: python
        Query().retrieve('[CH3]-[CH2X4]', ['BB-50-22Q1'])
    """

    enamine_dbs = ['BB-ForSale-22Q1', 'MADE-BB-23Q1-770M', 'REAL-Database-22Q1']

    def __init__(self, base_url: str = 'https://arthor.docking.org/'):
        self.base_url = base_url

    @property
    def dbs(self):
        return pd.DataFrame(requests.get(self.base_url + 'dt/data').json())

    def retrieve(self, query: str, dbnames: List[str]):
        """
        Returns a dataframe of the results of the query,
        with fields:

        * N_RB: number of rotatable bonds
        * N_HA: number of heavy atoms

        :param query: SMARTS query
        :param dbnames: list of names (see self.dbs)
        :return:
        """
        dbname: str = ','.join(dbnames)
        if isinstance(query, Chem.Mol):
            query = Chem.MolToSmarts(query)
        response: requests.Response = requests.get(self.base_url + f'/dt/{dbname}/search',
                                                   dict(query=query,
                                                        type='SMARTS',
                                                        length=1_000_000)
                                                   )
        data: dict = response.json()
        if data.get("message", '') == "SMARTS query is always false!":
            warnings.warn(f"SMARTS query {query} is always false")
            return pd.DataFrame()
        if data.get('warning',''):
            warnings.warn(data['warning'])
        if not data.get('recordsTotal', False):
            warnings.warn(f"SMARTS query {query} returned no matches")
            return pd.DataFrame()
        matches = pd.DataFrame(data['data'], columns=['idx', 'smiles_id', 'empty', 'something', 'db'])
        if len(matches) == 0:  # empty
            return matches
        matches['id'] = matches.smiles_id.str.split(expand=True)[1]
        matches['smiles'] = matches.smiles_id.str.split(expand=True)[0]
        matches.drop_duplicates('id')
        PandasTools.AddMoleculeColumnToFrame(matches, 'smiles', 'mol', includeFingerprints=True)
        matches = matches.loc[~matches.mol.isnull()]
        # tabs?
        matches['db'] = matches.db.str.strip()
        matches['N_RB'] = matches.mol.apply(AllChem.CalcNumRotatableBonds)
        matches['N_HA'] = matches.mol.apply(AllChem.CalcNumHeavyAtoms)
        return matches.sort_values('N_HA').reset_index(drop=True)
