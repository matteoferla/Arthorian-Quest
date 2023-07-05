from rdkit import Chem
from typing import Mapping, Sequence


def querimonate(mol: Chem.Mol,
                replacements: Mapping[int, str] = {},
                rgroups: Sequence[int] = {},
                generic_arocarbons: bool = False) -> Chem.Mol:
    """
    Given a molecule, convert it to a molecule with query atoms,
    with correct element, number of hydrogens and charge,
    but with overridden values as given by ``replacements`` argument,
    which accepts a dictionary of index (int) to SMARTS (str).

    A ``Chem.QueryAtom`` is a special atom with encoded ambiguity.
    ``Chem.MolFromSmarts`` will create a molecule with query atoms, but not a regular mol.
    A query atom has the additional methods
    ``.HasQuery()``, ``.GetQueryType()``, ``.DescribeQuery()``.
    cannot be instantiated from Python, but can be using ``Chem.AtomFromSmarts``.

    Additionally, any atom idx in argument ``rgroups`` will get a sequential isotope number from one
    and property 'R-group' of R + off-by-one number.
    If this index was not in replacements, then the SMARTS will have one connection more than there are
    and one implict hydrogen less.

    This function requires ``mol`` to have implicit hydrogens.

    ..code-block::python
       queried:Chem.Mol = querimonate(Chem.MolFromSmiles('c1cnccc1'), {2: '[c,n]'})
       Chem.MolToSmarts(queried)
       # '[c&H1]1:[c&H1]:[c,n]:[c&H1]:[c&H1]:[c&H1]:1'

    Note 1. ``atom.GetSmarts()`` is a method, but it could return '[N+]' already,
    which is complicated to deal with as appending at given positions may muck things up.
    And certain SMILES are not legal, for example 'CC[O+H2]' should be 'CC[OH2+]'

    Note 2. There is no OED word for 'to make into a query'... querify? Shmeh.
    Querimony is a complaint. There is no verb form. This is a jocular neologism,
    as RDKit will complain...
    """
    mod = Chem.RWMol(mol)
    atom: Chem.Atom
    for atom in mod.GetAtoms():
        assert atom.GetNumRadicalElectrons() == 0, 'This molecule has a radical'
        # if not isinstance(atom, Chem.QueryAtom):
        idx: int = atom.GetIdx()
        n_Hs: int = atom.GetImplicitValence() + atom.GetNumExplicitHs()
        symbol: str = atom.GetSymbol().lower() if atom.GetIsAromatic() else atom.GetSymbol()
        scharge: str = get_charge_string(atom)
        # pick relevant SMARTS
        if idx in replacements:
            smarts: str = replacements[idx]
        elif idx in rgroups:
            assert n_Hs == 0, 'R-group requested for zero Hs. Use charge or change element via ``replacement``'
            n_Xs: int = len(atom.GetNeighbors())
            smarts = f'[{symbol}H{n_Hs - 1}X{n_Xs + 1}]'
        elif generic_arocarbons and symbol == 'c' and scharge == '':
            smarts = f'[cH{n_Hs},nH0,oH0,sH0]'
        else:
            smarts = f'[{symbol}H{n_Hs}{scharge}]'
        # swap
        if isinstance(atom, Chem.QueryAtom):
            # weird...
            atom.SetQuery(smarts)
        else:
            mod.ReplaceAtom(idx, Chem.AtomFromSmarts(smarts), preserveProps=True)
    for r, idx in enumerate(rgroups):
        atom = mod.GetAtomWithIdx(idx)
        atom.SetIsotope(r + 1)
        atom.SetProp('R-group', f'R{r + 1}')
    return mod.GetMol()


def get_charge_string(atom: Chem.Atom) -> str:
    # charge
    charge: int = atom.GetFormalCharge()
    if charge > 0:
        return '+' * abs(charge)
    elif charge < 0:
        return '-' * abs(charge)
    else:
        return ''
