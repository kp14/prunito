import datetime
import io
import re
import pytest
from prunito.uniprot import (map_to_or_from_uniprot, search, current_release,
                             search_reviewed, search_unreviewed)
from prunito.utils import NoDataError, ExcessiveDataError


def test_map_to_or_from_uniprot_invalid_format():
    with pytest.raises(ValueError):
        r = map_to_or_from_uniprot(['2KLE'], 'pdb', 'acc') # has to be pdb_id


def test_map_to_or_from_uniprot_no_uniprot_specified():
    with pytest.raises(ValueError) as e_info:
        r = map_to_or_from_uniprot(['2KLE'], 'pdb_id', 'embl')
    assert e_info.value.args[0] == 'One of source or target format has to be UniProt ACC or ID.'


def test_map_to_or_from_uniprot():
    r = map_to_or_from_uniprot(['2KLE'], 'pdb_id', 'acc')
    assert r.map['2KLE'] == ['Q13563']


def test_current_release():
    regex = '20[0-9]{2}_[0,1][1-9]'
    r = current_release()
    assert re.match(regex, r)

def test_current_release_as_attribute():
    r = current_release()
    assert len(r.release()) == 7
    assert r.release()[4] == '_'
    assert r.release().startswith('20')


def test_current_release_as_date():
    r = current_release()
    assert isinstance(r.date(), datetime.datetime)


def test_search_reviewed():
    r = search_reviewed('name:tax-binding')
    assert r.size() == 1


# def test_search_reviewed_stringIO():
#     r = search_reviewed('name:tax-binding')
#     assert isinstance(r.fobject(), io.StringIO)


def test_search_unreviewed_no_reviewed_specified():
    r = search_unreviewed('taxonomy:191813 AND name:cytochrome')
    assert 'Unreviewed;' in r.text


def test_search_over_limit():
    with pytest.raises(ExcessiveDataError):
        _ = search('name:kinase')


def test_search_no_results():
    with pytest.raises(NoDataError):
        _ = search('name:kiniananase')


def test_search_wrong_format():
    with pytest.raises(ValueError):
        _ = search('name:tax', frmt='some')