from rdkit import Chem
from rdkit.Chem import AllChem


def test_embed_ethanol():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert AllChem.EmbedMolecule(mol, randomSeed=7) == 0
