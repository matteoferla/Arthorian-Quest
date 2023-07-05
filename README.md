# Arthorian Quest

> :construction: This is a work in progress. :construction:
> It will most likely be abandoned...

An experiment in using Arthor (arthor.docking.org) by NextMove and filtering the results with Fragmenstein,
in order to enumerate linkers of two hits, via a SMARTS pattern.

The key principle is that the SMARTS pattern can be a disconnected compound, thus allowing linkers to be searched for.
By creating a SMARTS pattern with a wished for expansion vector, one can search for linkers that connect two hits.

To wittle out impossible compounds, the results are filtered by 3D constrained conformer generation via Fragmenstein,
but using the SMARTS pattern making it lightning fast.

Step:

1. Generate a SMARTS pattern that reflects one's wishes for a molecule
2. Enumeration of possible compounds via Arthor
3. Filtering by 3D constrained conformer generation via Fragmenstein
4. Placement via Fragmenstein

## Generate SMARTS

In a SMARTS pattern an atom can be an element like a SMILES with explicit hydrogens,
say `[CH2]`, but can also be a comma separated list of alt elements say `[CH2,NH2]` (or better `[C,N;H2]`,
and number of interactions can can be specified with "X", say `[CH2X4]` means the carbon has 2 hydrogens and 4 connections.
If these are not specified, you have got yourself a connection vector, eh.
NB. In Arthor, the lowercase and uppercase R notations are ignored.

The function `querimonate` in [create.py](arthorian_quest/create.py) is intended to convert a regular `Chem.Mol`
to one that encodes SMARTS patterns (i.e. by having not `Chem.Atom` but `Chem.QueryAtom` instances).

```python
smol: Chem.Mol = querimonate(Chem.MolFromSmiles('c1cnccc1'), {2: '[c,n]'})
Chem.MolToSmarts(smol)
```
It's doc-string has this:

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

```python
queried:Chem.Mol = querimonate(Chem.MolFromSmiles('c1cnccc1'), {2: '[c,n]'})
Chem.MolToSmarts(queried)
# '[c&H1]1:[c&H1]:[c,n]:[c&H1]:[c&H1]:[c&H1]:1'
```

Note 1. ``atom.GetSmarts()`` is a method, but it could return '[N+]' already,
which is complicated to deal with as appending at given positions may muck things up.
And certain SMILES are not legal, for example 'CC[O+H2]' should be 'CC[OH2+]'

## Find hits

Querying [Arthor](arthor.docking.org) is simple and
is documented in the [Arthor API page](https://arthor.docking.org/api).

Nothing interesting here code-wise,
except the trick that the SMARTS can be a disconnected compound,
say `'C[OH1].NC(=O)N'`,
which is the basis for this whole project.

```python
from arthorian_quest import QueryArthor
from IPython.display import display
import pandas as pd

# ---- Available databases ----
display(QueryArthor().dbs)

# ---- Get query results ----
query='[aH0X3]([#6,#7,#8])1[c,n](-[Br,Cl,NH2])[aH0X2][aH0X2]a1.N-S(=O)(=O)-C'
results: pd.DataFrame = QueryArthor().retrieve(query=query, QueryArthor.enamine_dbs)
```

## Filtering by 3D

In [quick_place.py](arthorian_quest/quick_place.py) is `SmartMonsterHandler`,
which is not a subclass of Monster. As I mean to make it more in line with `Victor`.

A new monster object is created for each molecule because otherwise the parameters linker and cause interference.

There is something annoying that happens. The SMARTS pattern generated mapping may be poisonous, so gets ignored.
This is bad. As a result a check happens to prevent this. `monster.convert_origins_to_custom_map()` will have neg indices,
i.e. ignore this atom. When a molecule is ignored it's values are marked as failed.

```python
results = matches.mol.apply(SmartMonsterHandler([x1594, fippedSulfonamide], query, joining_cutoff=10))
matches['success'] = results.loc[~results.isna()].apply(operator.itemgetter('success'))
matches['ddG'] = results.loc[~results.isna()].apply(operator.itemgetter('ddG'))
matches['minimized_mol'] = results.loc[~results.isna()].apply(operator.itemgetter('mol'))
```

## Fragmenstein

The above needs to be incorporated in Victor, in order to get in protein results.
But actually for now. The results are rather surprising as the linkages are generally impossible/contorted,
regardless of protein.

## Etc.

## To do

* What happens when one submits a query that is always false?
* Make a better SMARTS generator
* Checks...
