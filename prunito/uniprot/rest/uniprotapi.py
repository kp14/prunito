import io
import requests
from collections import defaultdict

from ...utils import (UNIPROT_KNOWLEDGEBASE,
                               UNIPROT_BATCH,
                               UNIPROT_CONVERT,
                               UNIPROT_MAP,
                               VALID_ID_MAPPINGS,
                               is_value_in_iterable)

def current_release():
    """Get current public release of UniProtKB.

    Releases are specified as <year>_<release>, e.g. 2016_02.
    There are 11 or 12 releases per calendar year.

       Returns:
           str: Release, e.g. 2016_02

    """
    # Release number contained in request header
    # So we retrieve a random human entry and look up the value
    payload = {"query": "organism:9606 AND reviewed:yes",
               "random": "yes",
               }
    result = requests.get(UNIPROT_KNOWLEDGEBASE, params=payload)
    if result.ok:
        return result.headers['x-uniprot-release'] # Returns string by default
    else:
        result.raise_for_status()

def search_reviewed(query, frmt='txt', file=False):
    '''Search reviewed UniProtKB (Swiss-Prot) entries only.

    Accepts standard UniProtKB query syntax, so queries can be tested
    on the uniprot.org website.

    Args:
        query (str): UniProtKB query string.
        frmt (str; optional): Format for results.
            Can be txt, xml, rdf, fasta. Defaults to txt.
        file (bool): Whether to wrap returned string in StringIO.
            Defaults to False.

    Returns:
        str or None: Data, if any.
    '''
    result = _search(query, frmt=frmt, reviewed=True,
                     unreviewed=False, file=file)
    return result

def search_unreviewed(query, frmt='txt', file=False):
    '''Search unreviewed UniProtKB (TrEMBL) entries only.

    Accepts standard UniProtKB query syntax, so queries can be tested
    on the uniprot.org website.

    Args:
        query (str): UniProtKB query string.
        frmt (str; optional): Format for results.
            Can be txt, xml, rdf, fasta. Defaults to txt.
        file (bool): Whether to wrap returned string in StringIO.
            Defaults to False.

    Returns:
        str or None: Data, if any.
    '''
    result = _search(query, frmt=frmt, reviewed=False,
                     unreviewed=True, file=file)
    return result

def search_all(query, frmt='txt', file=False):
    '''Search all of UniProtKB (Swiss-Prot + TrEMBL).

    Accepts standard UniProtKB query syntax, so queries can be tested
    on the uniprot.org website.

    Args:
        query (str): UniProtKB query string.
        frmt (str; optional): Format for results.
            Can be txt, xml, rdf, fasta. Defaults to txt.
        file (bool): Whether to wrap returned string in StringIO.
            Defaults to False.

    Returns:
        str or None: Data, if any.
    '''
    result = _search(query, frmt=frmt, reviewed=True,
                     unreviewed=True, file=file)
    return result

def number_SP_hits(query):
    '''Search reviewed UniProtKB entries only (Swiss-Prot).

    Accepts standard UniProtKB query syntax.

    Args:
        query (str): UniProtKB query string.

    Returns:
        int: number of hits.
    '''
    result = _search(query, frmt='list', reviewed=True,
                     unreviewed=False, file=False)
    if result:
        hit_list = result.split('\n')
        number = len(hit_list) - 1
    else:
        number = 0
    return number


def retrieve_batch(ac_list, frmt='txt', file=False):
    '''Batch retrieval of uniProtKB entries.

    Args:
        ac_list (list): UniProtKB accessions
        frmt (str; optional): Format for results.
            Can be txt, xml, rdf, fasta. Defaults to txt.
        file (bool): Whether to wrap returned string in StringIO.
            Defaults to False.

    Returns:
        str or None: Data, if any.
    '''
    payload = {'query':' '.join(ac_list),
               'format':frmt}
    result = requests.get(UNIPROT_BATCH, params=payload)
    if result.ok:
        if len(result.content) > 0:
            if file:
                return _to_StringIO(result.content)
            else:
                return str(result.content, encoding="ascii")
        else:
            return None
    else:
        result.raise_for_status()


def convert(path, source_fmt='txt', target_fmt='xml', encoding='ascii'):
    '''Convert UniProt entries between different formats.

    Currently, only single entries are supported.

    Parameters:
        path (str): Path to file containing entry to be converted.
        source_fmt (str): Source format. Default: txt.
        target_fmt (str): Target format. Default: xml.

    Returns:
        str: Converted data.
    '''
    payload = {'type': 'uniprot',
               'from': source_fmt,
               'to': target_fmt,
               }
    files = {'data': open(path, 'r', encoding=encoding)}
    response = requests.post(UNIPROT_CONVERT,
                             data=payload,
                             files=files
                             )
    if response.ok:
        return response.text
    else:
        response.raise_for_status()


def map_to_or_from_uniprot(id_list, source_fmt, target_fmt):
    '''Map one set of identifiers to another.

    UniProt provides mappings to more than 100 databases it cross-references.
    The limitation of the service is that it only maps X -> UniProt or
    UniProt -> Y. A direct mapping from X -> Y is not possible; this would be
    a 2-step process.

    Mappings are not necessarily 1-to-1 as e.g. a given UniProtKB entry can
    have many associated PDB cross-references. To reflect this, the returned
    data structure is a dict.

    As mapped data are returned in tabular format, customization of columns is
    possible in principle. However, this is not implemented here as the function
    is only concerned with the actual mapping and not any further data. Those
    could be retrieved in a second call using search or batch retrieval.

    See http://www.uniprot.org/help/programmatic_access#conversion for details.
    Note: The response.url field contains the URL from which to download
    the mapping, e.g. http://www.uniprot.org/mapping/M20160504763V34ZKX0.tab

    Args:
        id_list (list): Identifiers to be mapped. Identifiers should be strings.
        source_fmt (str): Format of the provided identifiers. See UniProt help for allowed values.
        target_fmt (str): Desired format of the identifiers. See UniProt help for allowed values.

    Returns:
        dict: With source IDs as keys and a list of target IDs a values.

    Raises:
        ValueError: If invalid identifier formats are used.
        ValueError: If not at least source or target is UniProt accession/ID
    '''
    if source_fmt.upper() not in VALID_ID_MAPPINGS or target_fmt.upper() not in VALID_ID_MAPPINGS:
        msg = 'Invalid mapping source(s): {0}, {1}\nUse one of:\n{2}'.format(source_fmt.upper(),
                                                                             target_fmt.upper(),
                                                                             ', '.join(VALID_ID_MAPPINGS))
        raise ValueError(msg)
    if not source_fmt.upper() in ('ACC', 'ID'):
        if not target_fmt.upper() in ('ACC', 'ID'):
            raise ValueError('Source or target format has to be UniProt ACC or ID.')
    id_list = ' '.join(id_list)
    payload = {'from': source_fmt.upper(),
               'to': target_fmt.upper(),
               'format': 'tab',
               'query': id_list,
               }
    response = requests.get(UNIPROT_MAP, params=payload)
    if response.ok:
        mapped = defaultdict(list)
        lines = iter(response.text.split('\n'))
        lines.__next__() # skip header
        for line in lines:
            try:
                source, target = line.strip().split('\t')
            except ValueError:
                pass
            else:
                mapped[source].append(target)
        return mapped
    else:
        response.raise_for_status()


def _search(query, frmt='txt', reviewed=True, unreviewed=True, file=False):
    _check_format(frmt)
    payload = {'query':query, 'format':frmt}
    if reviewed and unreviewed:
        pass
    elif reviewed and not unreviewed:#Swiss-Prot
        payload['query'] += ' AND reviewed:yes'
    elif not reviewed and unreviewed:#TrEMBL
        payload['query'] += ' AND reviewed:no'
    elif not reviewed and not unreviewed:
        msg = ('At least one of parameters `reviewed` and `unreviewed` has to be True.\n'
               'Found: reviewed: {0}, unreviewed: {1}')
        raise ValueError(msg.format(reviewed, unreviewed))
    result = requests.get(UNIPROT_KNOWLEDGEBASE, params=payload)
    if result.ok:
        if len(result.content) > 0:
            if file:
                return _to_StringIO(result.content)
            else:
                return str(result.content, encoding="ascii")
        else:
            return None
    else:
        result.raise_for_status()


def _to_StringIO(text):
    return io.StringIO(text.decode())
    #return io.StringIO(unicode(text))

def _check_format(fmt):
    return_formats = ('html',
                      'tab',
                      'xls',
                      'fasta',
                      'gff',
                      'txt',
                      'xml',
                      'rdf',
                      'list',
                      #'rss',
                      )
    if not is_value_in_iterable(fmt, return_formats):
        msg = 'Allowed values: {0}\nPassed in value: {1}'
        raise ValueError(msg.format(return_formats, fmt))
