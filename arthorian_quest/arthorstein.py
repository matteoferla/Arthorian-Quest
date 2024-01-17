__all__ = ['prep_for_laboratory', 'run_experiment', 'create_laboratory_df', 'assess_experiment', 'get_custom_map']

from rdkit import Chem
from rdkit.Chem import AllChem
import contextlib
import pandas as pd
import warnings
from .query import QueryArthor
from typing import Dict, Union, Optional


def prep_for_laboratory(query_mol: Chem.Mol,
                        template_mol: Chem.Mol,
                        experiment_name:str='') \
        -> pd.DataFrame:
    """
    This function is a wrapper for ``create_laboratory_df`` and ``assess_experiment``.
    It creates a dataframe of analogues and assesses the experiment.
    The dataframe can be saved and run in a different node.
    The logic being that it is often required to tweak the SMARTS query to get a good amount of analogues,
    whereas the placement can happen in the background.
    To run one to check run ``assess_experiment``.

    :param query_mol:
    :param template_mol:
    :param experiment_name:
    :return:
    """

    analogs = create_laboratory_df(query_mol, template_mol, experiment_name)
    print('No analogues found' if len(analogs) == 0 else f'{len(analogs)} analogues found')
    return analogs
def run_experiment(query: Chem.Mol, template: Chem.Mol, pdb_block, experiment_name: str=''):
    """
    Calls ``create_laboratory_df`` and ``Laboratory.place``.

    :param query:
    :param template:
    :param pdb_block:
    :param experiment_name:
    :return:
    """
    from fragmenstein import Laboratory

    placements = pd.DataFrame()
    with contextlib.suppress(KeyError):
        tobeplaced = create_laboratory_df(query, template, experiment_name=experiment_name)
        placements = Laboratory(pdb_block).place(tobeplaced)
    placements['experiment'] = experiment_name
    return placements

def create_laboratory_df(query_mol: Chem.Mol, template_mol: Chem.Mol, experiment_name: str= '') -> pd.DataFrame:
    """
    Given a query molecule and a template molecule, return a dataframe of analogues
    by searching Arthor.
    This dataframe is ready to be passed to ``Laboratory.place``.
    It has the columns:

    * name: the name of the analogue
    * hits: a list of the template molecule
    * custom_map: a dictionary of the custom map from the template to the analogue
    * experiment: the name of the experiment

    :param query_mol:
    :param template_mol:
    :param experiment_name:
    :return:
    """
    # ## Sanitize
    if not experiment_name and query_mol.HasProp('experiment'):
        experiment_name = query_mol.GetProp('experiment')
    query_mol.SetProp('experiment', experiment_name)
    query_mol.SetProp('template', template_mol.GetProp('_Name'))
    # ## Run search
    arthor = QueryArthor()
    df = arthor.retrieve(Chem.MolToSmarts(query_mol), ['real-database-22q1'])
    query_mol.UpdatePropertyCache()
    #AllChem.SanitizeMol(query)
    # exit early if no hits
    if len(df) == 0:
        return pd.DataFrame()
    try:
        df = df.rename(columns={'id': 'name'})
        df['hits'] = df.smiles.apply(lambda s: [template_mol])
        mapper = lambda smiles: {template_mol.GetProp('_Name'): get_custom_map(query_mol, template_mol, smiles)}
        df['custom_map'] = df.smiles.apply(mapper)
        df['experiment'] = experiment_name
        return df
    except Exception as error:
        print(error.__class__.__name__, error)
        return df

def assess_experiment(query: Chem.Mol, template: Chem.Mol, pdb_block: str) -> Chem.Mol:
    """
    Given a query molecule and a template molecule, assess the experiment by placing the query molecule

    :param query:
    :param template:
    :param pdb_block:
    :return:
    """
    from fragmenstein import Laboratory

    name = query.GetProp('experiment') if query.HasProp('experiment') else 'test'
    query.UpdatePropertyCache()
    AllChem.SanitizeMol(query, catchErrors=True)
    assert query.GetConformer()
    cm = get_custom_map(query, template, Chem.MolToSmiles(template))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        minimised = Laboratory.Victor([template], pdb_block=pdb_block).place(Chem.MolToSmiles(template),
                                                                        long_name=name,
                                                                        custom_map=cm).minimized_mol
        assert minimised, 'Nothing returned'
        return minimised



def get_custom_map(shared_query: Chem.Mol, template: Chem.Mol, final: Union[str, Chem.Mol]) -> Dict[int, int]:
    """
    Given a template molecule and a final molecule and a shared query molecule,
    return the custom map going from atom indices in the template molecule to the final molecule in the smiles string.

    :param shared_query:
    :param template:
    :param final:
    :return:
    """
    query2template = dict(enumerate(template.GetSubstructMatch(shared_query)))
    if isinstance(final, str):
        final = Chem.MolFromSmiles(final)
    query2analog_iter = enumerate(final.GetSubstructMatch(shared_query))
    return {query2template[qi]: ai for qi, ai in query2analog_iter}
