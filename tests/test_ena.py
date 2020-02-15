import pytest
import requests
from prunito import ena
from prunito.utils import NoDataError


def test_retrieve_fasta():
    r = ena.retrieve("BAC67592", fmt="fasta")
    assert r.startswith(">ENA|BAC67592")


def test_retrieve_xml():
    r = ena.retrieve("BAC67592", fmt="xml")
    assert r.startswith('<?xml version="1.0" encoding="UTF-8"?>')


def test_retrieve_fasta_wrong_id():
    with pytest.raises(NoDataError):
        r = ena.retrieve("BAC00000", fmt="fasta")


def test_retrieve_embl():
    r = ena.retrieve("BAC67592", fmt="text")
    assert r.startswith("ID   BAC67592;")


def test_retrieve_invalid_id():
    with pytest.raises(ValueError):
        r = ena.retrieve("111BAC67592")


def test_retrieve_many_ids():
    with pytest.raises(ValueError):
        r = ena.retrieve("BAC67592 BAC67593 BAC67580")


def test_retrieve_wrong_url():
    with pytest.raises(ValueError):
        r = ena.retrieve("display=111BAC6&7592=some")


def test_retrieve_wrong_url_xml():
    with pytest.raises(ValueError):
        r = ena.retrieve("display=111BAC6&7592=some", fmt="xml")


def test_translate_start_job(monkeypatch):
    pass
