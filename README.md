# arthorian-quest

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

The 

```python
results = matches.mol.apply(SmartMonsterHandler([x1594, fippedSulfonamide], query, joining_cutoff=10))
matches['success'] = results.loc[~results.isna()].apply(operator.itemgetter('success'))
matches['ddG'] = results.loc[~results.isna()].apply(operator.itemgetter('ddG'))
matches['minimized_mol'] = results.loc[~results.isna()].apply(operator.itemgetter('mol'))
```


## Etc.

## To do

* What happens when one submits a query that is always false?
* Checks...
