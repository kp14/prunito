# -*- coding: utf-8 -*-
"""Module for accessing ENA data via REST-ful web services."""
import re
import requests
from .utils import WSResponse


from .utils import ENA_DATA, NoDataError, ENA_IDENTIFIER_REGEX

session = requests.Session()


def retrieve(identifier, fmt='fasta'):
    '''Retrieve data based on an ENA identifier.

    This is intended for ENA accession numbers and protein IDs but I guess it
    should work for assembly metadata (e.g. GCA_001521735) or BioProject metadata
    (PRJNA301708), too. Not sure how useful the latter are though.

    Args:
        identifier (str): Identifier to retrieve.
        fmt (str, optional): Retrieval format. Can be fasta, text, xml. Defaults to fasta.

    Returns:
        WSResponse

    Raises:
        NoDataError: If no data are returned
    '''
    match = re.match(ENA_IDENTIFIER_REGEX, identifier)
    if not match or len(match.group()) is not len(identifier):
        raise ValueError('Identifier seems to have wrong format. Provide only 1 identifier.')
    error_msg = 'display type is either not supported or entry is not found'
    data = identifier + '&display=' + fmt
    result = WSResponse(session.get('/'.join([ENA_DATA, data])))
    if result.ok:
        if error_msg in result.text:
            raise NoDataError(error_msg)
        return result
    else:
        result.raise_for_status()
