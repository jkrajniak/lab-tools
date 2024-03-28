from src.md_tools import monomer_purger


def test_clean_atoms():
    atoms = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    monomers = [2, 4, 6, 8]
    mapping = monomer_purger.clean_atoms(atoms, monomers)
    assert mapping == {1: 1, 3: 2, 5: 3, 7: 4, 9: 5}