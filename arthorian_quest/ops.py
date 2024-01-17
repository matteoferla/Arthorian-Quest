__all__ = ['silence_rdkit', 'show_experiment', 'retrieve_smartsplus']

from rdkit import RDLogger, Chem
from rdkit.Chem import AllChem, Draw
from IPython.display import display
from typing import Union
import PIL
import requests
import io


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



def retrieve_smartsplus(smarts: Union[str, Chem.Mol], PIL_image=True, **options) -> Union[display.Image, PIL.Image]:
    """
    Given a SMARTS query, retrieve the image from https://smarts.plus.
    The returned object is an IPython.display.Image not a PIL.Image.
    If using this image remember to cite it as
    "SMARTSviewer smartsview.zbh.uni-hamburg.de, ZBH Center for Bioinformatics, University of Hamburg"

    :param smarts: SMARTS query or Chem.Mol
    :param PIL_image: return PIL.Image instead of IPython.display.Image
    :param options: See https://smarts.plus/rest
    :return:
    """
    if isinstance(smarts, Chem.Mol):
        q = smarts
        smarts: str = Chem.MolFromSmarts(q)
    # retrieve from smarts.plus
    response: requests.Response = requests.get('https://smarts.plus/smartsview/download_rest', # noqa to sensitive data in options
                                               {'smarts': smarts, **options}
                                              )
    png_binary: bytes = response.content
    if PIL_image:
        return PIL.Image.open(io.BytesIO(png_binary))
    else:
        return display.Image(data=png_binary)

