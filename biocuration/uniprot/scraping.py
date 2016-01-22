# -*- coding: utf-8 -*-
"""
Some UniProt data are not available via the REST API, these can be scraped.

Created on Fri Jan 22 13:04:00 2016

@author: kp14
"""
import requests
from bs4 import BeautifulSoup

from ..utils import UNIPROT_TAXONOMY


session = requests.Session()


def get_lineage(taxid, full=True):
    '''Retrieve taxonomic lineage of given taxid.

    Parameter:
        taxid: taxid string, e.g. `9606`
        full: boolean; UniProtKB entries contain a shortened lineage which can
                       be retrieved by setting full=False
    Returns:
        list of taxon nodes (names)
    '''
    page = session.get('/'.join([UNIPROT_TAXONOMY, taxid]))
    if page.ok:
        decoded = page.content.decode('utf8')
        soup = BeautifulSoup(decoded, 'html.parser')
        # Using this attribute I figuered out by experimenting
        arefs = soup.findAll(property='rdfs:subClassOf')
        # text: contains the node name, attrs[ref] the taxid
        if full:
            return [(ref.text, ref.attrs['href']) for ref in arefs]
        return [(ref.text, ref.attrs['href']) for ref in arefs if not 'style' in ref.attrs]
