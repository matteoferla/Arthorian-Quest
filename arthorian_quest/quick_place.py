__all__ = ['SmartMonsterHandler']

from rdkit import Chem
from fragmenstein import Monster
from contextlib import suppress
from rdkit.rdBase import BlockLogs


class SmartMonsterHandler:
    """
    Given hits and a SMARTS query w/ dot separated patterns,
    wherein each sequentially matches the hits

    ... code-block:: python
        sm = SmartMonsterHandler([hit1, hit2], query, joining_cutoff=10)
        sm(matches.mol[0])

    or

    ... code-block:: python
        sm = SmartMonsterHandler([hit1, hit2], query, joining_cutoff=10)
        results = matches.mol.apply(sm)
        matches['success'] = results.loc[~results.isna()].apply(operator.itemgetter('success'))
        matches['ddG'] = results.loc[~results.isna()].apply(operator.itemgetter('ddG'))
        matches['minimized_mol'] = results.loc[~results.isna()].apply(operator.itemgetter('mol'))

"""

    def __init__(self, hits, query, **monster_options):
        self.hits = hits
        self.monster_options = monster_options
        self.names = [h.GetProp('_Name') for h in self.hits]
        self.smols = [Chem.MolFromSmarts(q) for q in query.split('.')]
        self.mappings = [h.GetSubstructMatch(smol) for h, smol in zip(self.hits, self.smols)]
        assert all([len(m) > 0 for m in self.mappings]), 'The query does not match the hits'
        self._exception = Exception

    def map(self, mol):
        return {n: dict(zip(m, mol.GetSubstructMatch(s))) for n, m, s in zip(self.names, self.mappings, self.smols)}

    def __call__(self, mol):
        monster = Monster(self.hits, **self.monster_options)
        mapping = self.map(mol)
        with suppress(self._exception), BlockLogs():
            monster.place(mol, custom_map=mapping)
            used_map = monster.convert_origins_to_custom_map()
            success = monster.mmff_minimize(allow_lax=False)
            ddG = monster.MMFF_score()
            if any([len(set(mapping[name].items()) - set(used_map[name].items())) != 0 for name in self.names]):
                success = False
                ddG = float('nan')
            return {'mol': monster.positioned_mol,
                    'custom_map': mapping,
                    'success': success,
                    'ddG': ddG
                    }
        return {'mol': mol, 'custom_map': mapping, 'success': False, 'ddG': float('nan')}