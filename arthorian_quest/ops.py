from .query import QueryArthor

import contextlib
import pandas as pd
import warnings
from rdkit import RDLogger, Chem
from rdkit.Chem import AllChem, Draw
from IPython.display import display
from fragmenstein import Laboratory, Victor, Wictor

logger = RDLogger.logger()
logger.setLevel(RDLogger.CRITICAL)
Laboratory.Victor = Wictor

to_run = {}

def get_custom_map(query: Chem.Mol, template: Chem.Mol, smiles: str):
    name = template.GetProp('_Name')
    query2template = dict(enumerate(template.GetSubstructMatch(query)))
    analog = Chem.MolFromSmiles(smiles)
    query2analog_iter = enumerate(analog.GetSubstructMatch(query))
    return {name: {query2template[qi]: ai for qi, ai in query2analog_iter}}

def prep_for_place(query: Chem.Mol, template: Chem.Mol, experiment_name: str='') -> pd.DataFrame:
    arthor = QueryArthor()
    df = arthor.retrieve(Chem.MolToSmarts(query), ['real-database-22q1'])
    query.UpdatePropertyCache()
    #AllChem.SanitizeMol(query)
    if not experiment_name:
        experiment_name = query.GetProp('experiment') if query.HasProp('experiment') else 'test'
    if len(df) == 0:
        return pd.DataFrame()
    try:
        df = df.rename(columns={'id': 'name'})
        df['hits'] = df.smiles.apply(lambda s: [template])
        df['custom_map'] = df.smiles.apply(lambda s: get_custom_map(query, template, s))
        df['experiment'] = experiment_name
        return df
    except Exception as error:
        print(error.__class__.__name__, error)
        return df

def show_experiment(query3d: Chem.Mol, experiment_name: str=''):
    query = Chem.Mol(query3d)
    if not experiment_name:
        experiment_name = query.GetProp('experiment') if query.HasProp('experiment') else 'test'
    AllChem.Compute2DCoords(query)
    print(Chem.MolToSmarts(query))
    display(query)
    Draw.MolToFile(query, f'{experiment_name}.png')
    return None

def run_experiment(query: Chem.Mol, template: Chem.Mol, pdb_block, experiment_name: str=''):
    if not experiment_name:
        experiment_name = query.GetProp('experiment') if query.HasProp('experiment') else 'test'
    placements = pd.DataFrame()
    with contextlib.suppress(KeyError):
        tobeplaced = prep_for_place(query, template, experiment_name)
        Laboratory(pdb_block).place(tobeplaced)
    placements['experiment'] = experiment_name
    print(f'{len(placements)} found')
    return placements

def assess_experiment(query: Chem.Mol, template: Chem.Mol, pdb_block: str):
    name = query.GetProp('experiment') if query.HasProp('experiment') else 'test'
    query.UpdatePropertyCache()
    AllChem.SanitizeMol(query, catchErrors=True)
    assert query.GetConformer()
    cm = get_custom_map(query, template, Chem.MolToSmiles(template))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert Wictor([template], pdb_block=pdb_block).place(Chem.MolToSmiles(template), long_name=name, custom_map=cm).minimized_mol

def boilerplace(query, template, experiment, pdb_block: str) -> pd.DataFrame:
    query.SetProp('experiment', experiment)
    query.SetProp('template', template.GetProp('_Name') )
    assess_experiment(query, template, pdb_block)
    show_experiment(query)
    analogs = prep_for_place(query, template, experiment)
    print('No analogues found' if len(analogs) == 0 else f'{len(analogs)} analogues found')
    return analogs
