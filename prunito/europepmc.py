"""Module for accessing EuropePMC via REST-ful web services."""

import requests
import urllib
from collections import OrderedDict
from .utils import EPMC_SEARCH


def search(query, fmt='json', result_type='lite', page_size=25):
    """Search publication database.

    The details of the search syntax are explain here:
    http://europepmc.org/Help#SSR

    The parameter pageSize is provided which allows downloading of at
    maximum 1000 hits, enough for our purposes. Thus, a parameter page
    or cursorMark are not deemed necessary.

    Args:
        query (str): The Query.
        fmt (str, optional): Format of returned data. Can be XML, JSON, DC.
            Defaults to JSON.
        result_type (str, optional): Depth of returned data.
            Can be idlist, core, lite.
            idlist: returns a list of IDs and sources for the
                given search terms
            lite: returns key metadata for the given search terms (default).
            core: returns full metadata for a given
                publication ID; including abstract, full text
                links, and MeSH terms.
        pageSize (int, optional): Specify the number of articles you wish to
            retrieve in each page. Max: 1000, Default: 25.
    """
    # Normal dict did not work as parameters were shuffled
    payload = OrderedDict()
    payload['query'] = query
    payload['format'] = fmt
    payload['resulttype'] = result_type
    payload['pageSize'] = page_size

    result = requests.get(EPMC_SEARCH + urllib.parse.urlencode(payload))
    if result.ok:
        if fmt == 'json':
            return result.json()
        else:
            return result.text


def get_pmid_metadata(pmid):
    """Retrieve core metadata of a given PubMed ID.

    Core metadata include the abstract, full text links and MeSH terms.
    All of these have to be extracted from the JSON.

    As a PubMed ID should be unique, this function returns zero or one
    result. If there are more than one hits then the first one is returned.

    Args:
        pmid (str or int): PubMed ID.

    Returns:
        Metadata (dict) or None
    """
    query = 'ext_id:{} src:med'.format(str(pmid))
    result = search(query, result_type='core')
    if result['hitCount'] == 0:
        return None
    else:
        return result['resultList']['result'][0]
