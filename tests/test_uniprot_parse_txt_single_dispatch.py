import os
import pathlib
import pytest
from io import StringIO
from prunito import uniprot as up


# Import paths for data files
base_dir = os.path.dirname(__file__)
filename = 'one_sp_entry.txl'
TEST_DATA_LOCATION = os.path.join(base_dir, 'SwissProt', filename)


EXPECTED_ACC = 'Q6E804'


def test_parser_with_file_object():
    with open(TEST_DATA_LOCATION, 'r') as infile:
        for record in up.parse_txt(infile):
            assert record.primary_accession == EXPECTED_ACC


def test_parser_with_StringIO():
    with open(TEST_DATA_LOCATION, 'r') as infile:
        content = infile.read()
        stringio = StringIO(content)
    for record in up.parse_txt(stringio):
        assert record.primary_accession == EXPECTED_ACC


def test_parser_with_Path_instance():
    p = pathlib.Path(TEST_DATA_LOCATION)
    for record in up.parse_txt(p):
        assert record.primary_accession == EXPECTED_ACC


def test_parser_with_str_path():
    for record in up.parse_txt(TEST_DATA_LOCATION):
        assert record.primary_accession == EXPECTED_ACC


def test_parser_with_unsupported_parameter_int():
    with pytest.raises(NotImplementedError):
        for record in up.parse_txt(1):
            pass


def test_parser_with_unsupported_parameter_dict():
    d = {1: 'some'}
    with pytest.raises(NotImplementedError):
        for record in up.parse_txt(d):
            pass