import pytest
from prunito.uniprot import map_to_or_from_uniprot


def test_map_to_or_from_uniprot_invalid_format():
    with pytest.raises(ValueError):
        r = map_to_or_from_uniprot(['2KLE'], 'pdb', 'acc') # has to be pdb_id


def test_map_to_or_from_uniprot_no_uniprot_specified():
    with pytest.raises(ValueError) as e_info:
        r = map_to_or_from_uniprot(['2KLE'], 'pdb_id', 'embl')
    assert e_info.value.args[0] == 'Source or target format has to be UniProt ACC or ID.'


def test_map_to_or_from_uniprot():
    r = map_to_or_from_uniprot(['2KLE'], 'pdb_id', 'acc')
    assert isinstance(r, dict)
    assert r['2KLE'] == ['Q13563']