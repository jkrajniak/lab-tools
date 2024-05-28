from md_tools.md_libs.files_io import LammpsReader


def test_get_pandas_dataframe():
    lmp = LammpsReader()
    lmp.read_data('tests/data.polymer')

    df = lmp.get_atoms_as_dataframe()
    assert df.shape == (100, 7)
    assert set(df.columns) == {'atom_type', 'res_id', 'position', 'image', 'charge', 'vel', 'mass'}
