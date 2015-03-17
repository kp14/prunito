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
import sys


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


def _search(query, fmt='json', resulttype='lite', page='1',
            dataset='metadata', callback=None):
    """Does all the heavy lifting.

    JSON is always returned.
    """
    base_url = "http://www.ebi.ac.uk/europepmc/webservices/rest/search/"
    payload = {"query": query,
               "format": fmt,
               "resulttype": resulttype,
               "page": page,
               "dataset": dataset}
    if callback is not None:
        payload['callback'] = callback
    params = _urlencode_no_plus(payload)
    # Cannot use requests package as:
    # [1] Requests prepends a '?' to the query string which europmc
    # does not like. [2] Even providing a preformed URL to requests
    # doesn't solve this as requests then calls again urlencode
    # which srews everything up. So we go for urllib.urlopen()
    url = '{0}{1}'.format(base_url, params)
    print url
    result = urllib.urlopen(url)
    data = result.read()
    return data


def _urlencode_no_plus(query, doseq=False, safe='', encoding=None, errors=None):
    """Encode a dict or sequence of two-element tuples into a URL query string.

    If any values in the query arg are sequences and doseq is true, each
    sequence element is converted to a separate parameter.

    If the query arg is a sequence of two-element tuples, the order of the
    parameters in the output will match the order of parameters in the
    input.

    The components of a query arg may each be either a string or a bytes type.
    When a component is a string, the safe, encoding and error parameters are
    sent to the urllib.quote function for encoding.
    """

    if hasattr(query, "items"):
        query = query.items()
    else:
        # It's a bother at times that strings and string-like objects are
        # sequences.
        try:
            # non-sequence items should not work with len()
            # non-empty strings will fail this
            if len(query) and not isinstance(query[0], tuple):
                raise TypeError
            # Zero-length sequences of all types will get here and succeed,
            # but that's a minor nit.  Since the original implementation
            # allowed empty dicts that type of behavior probably should be
            # preserved for consistency
        except TypeError:
            ty, va, tb = sys.exc_info()
            raise TypeError("not a valid non-string sequence "
                            "or mapping object").with_traceback(tb)

    l = []
    if not doseq:
        for k, v in query:
            if isinstance(k, bytes):
                k = urllib.quote(k, safe)
            else:
                k = urllib.quote(str(k), safe, encoding, errors)

            if isinstance(v, bytes):
                v = urllib.quote(v, safe)
            else:
                v = urllib.quote(str(v), safe, encoding, errors)
            l.append(k + '=' + v)
    else:
        for k, v in query:
            if isinstance(k, bytes):
                k = urllib.quote(k, safe)
            else:
                k = urllib.quote(str(k), safe, encoding, errors)

            if isinstance(v, bytes):
                v = urllib.quote(v, safe)
                l.append(k + '=' + v)
            elif isinstance(v, str):
                v = urllib.quote(v, safe, encoding, errors)
                l.append(k + '=' + v)
            else:
                try:
                    # Is this a sufficient test for sequence-ness?
                    x = len(v)
                except TypeError:
                    # not a sequence
                    v = urllib.quote(str(v), safe, encoding, errors)
                    l.append(k + '=' + v)
                else:
                    # loop over the sequence
                    for elt in v:
                        if isinstance(elt, bytes):
                            elt = urllib.quote(elt, safe)
                        else:
                            elt = urllib.quote(str(elt), safe, encoding, errors)
                        l.append(k + '=' + elt)
    return '&'.join(l)

# Utilities to parse URLs (most of these return None for missing parts):
# unwrap('<URL:type://host/path>') --> 'type://host/path'
# splittype('type:opaquestring') --> 'type', 'opaquestring'
# splithost('//host[:port]/path') --> 'host[:port]', '/path'
# splituser('user[:passwd]@host[:port]') --> 'user[:passwd]', 'host[:port]'
# splitpasswd('user:passwd') -> 'user', 'passwd'
# splitport('host:port') --> 'host', 'port'
# splitquery('/path?query') --> '/path', 'query'
# splittag('/path#tag') --> '/path', 'tag'
# splitattr('/path;attr1=value1;attr2=value2;...') ->
#   '/path', ['attr1=value1', 'attr2=value2', ...]
# splitvalue('attr=value') --> 'attr', 'value'
# urllib.parse.unquote('abc%20def') -> 'abc def'
# quote('abc def') -> 'abc%20def')
