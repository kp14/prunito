"""Module for accessing EuropePMC via REST-ful web services.

This uses the `search` module of said web services.
URL constructing follows the following pattern:

GET http://www.ebi.ac.uk/europepmc/webservices/rest/search/
                            query={}[&parameters]

Mandatory parameter:

`query`: various fields can be searched

Optional parameters:

`resulttype`: idlist, core, lite

            idlist: returns a list of IDs and sources for the
                given search terms
            lite: returns key metadata for the given search terms (default).
            core: returns full metadata for a given
                publication ID; including abstract, full text
                links, and MeSH terms.dataset: metadata,
                fulltext (default=metadata)

`page`: Specify the results page you wish to retrieve, where
        applicable. Page length of 25 and page numbers start
        at 1; default value is 1 (i.e. publications 1-25), if
        parameter unspecified.

`format`: XML, JSON (default=JSON)

`callback`: Use the callback parameter to make cross-domain
        requests and retrieve JSON data. The value is user
        specified. The web service wraps its response with
        the given callback value.
"""
import requests
import urllib

from collections import OrderedDict


def search(query, fmt='json', resulttype='lite', page='1',
           dataset='metadata', callback=None):
    """Search with default values for all parameters.
    Returns JSON data.
    """
    result = _search(query, fmt=fmt, resulttype=resulttype,
                     page=page, dataset=dataset, callback=callback)
    return result


def search_fulltext(query, fmt='json', resulttype='',
                    dataset='fulltext', page='1', callback=None):
    """Search fulltext articles only.
    Returns JSON data.
    """
    result = _search(query, fmt=fmt, resulttype=resulttype,
                     page=page, dataset=dataset, callback=callback)
    return result


def retrieve_pmid(pmid):
    query = 'ext_id:{}'.format(str(pmid))
    result = _search(query, fmt='json', resulttype='core',
                     page='1', dataset='metadata', callback=None)
    return result


def _search(query, fmt='json', resulttype='lite', page='1',
            dataset='metadata', callback=None):
    """Does all the heavy lifting.

    JSON is always returned.
    """
    base_url = "http://www.ebi.ac.uk/europepmc/webservices/rest/search/"

    # Normal dict did not work as parameters were shuffled
    payload = OrderedDict()
    payload['query'] = query
    payload['format'] = fmt
    payload['resulttype'] = resulttype
    payload['page'] = page
    payload['dataset'] = dataset

    if callback is not None:
        payload['callback'] = callback

    result = requests.get(base_url + urllib.parse.urlencode(payload))
    if result.ok:
        if fmt == 'json':
            return result.json()
        else:
            return result.text
