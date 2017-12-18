import io
import requests
from collections import defaultdict

from ...utils import (UNIPROT_KNOWLEDGEBASE,
                      UNIPROT_BATCH,
                      UNIPROT_CONVERT,
                      UNIPROT_MAP,
                      VALID_ID_MAPPINGS,
                      NoDataError,
                      )


session = requests.Session()


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
               'format': 'list',
               }
    result = session.get(UNIPROT_KNOWLEDGEBASE, params=payload)
    if result.ok:
        return result.headers['x-uniprot-release'] # Returns string by default
    else:
        result.raise_for_status()


def search_reviewed(query, frmt='txt', file_handle=False):
    '''Search reviewed UniProtKB (Swiss-Prot) entries only.

    Accepts standard UniProtKB query syntax, so queries can be tested
    on the uniprot.org website.

    Args:
        query (str): UniProtKB query string.
        frmt (str; optional): Format for results.
            Can be txt, xml, rdf, fasta. Defaults to txt.
        file_handle (bool): Whether to wrap returned string in StringIO.
            Defaults to False.

    Returns:
        str or None: Data, if any.
    '''
    if 'reviewed:yes' not in query:
        query += ' AND reviewed:yes'
    result = _search(query, frmt=frmt, file_handle=file_handle)
    return result


def search_unreviewed(query, frmt='txt', file_handle=False):
    '''Search unreviewed UniProtKB (TrEMBL) entries only.

    Accepts standard UniProtKB query syntax, so queries can be tested
    on the uniprot.org website.

    Args:
        query (str): UniProtKB query string.
        frmt (str; optional): Format for results.
            Can be txt, xml, rdf, fasta. Defaults to txt.
        file_handle (bool): Whether to wrap returned string in StringIO.
            Defaults to False.

    Returns:
        str or None: Data, if any.
    '''
    if 'reviewed:no' not in query:
        query += ' AND reviewed:no'
    result = _search(query, frmt=frmt, file_handle=file_handle)
    return result


def search(query, frmt='txt', file_handle=False):
    '''Search UniProtKB.

    Accepts standard UniProtKB query syntax, so queries can be tested
    on the uniprot.org website. Unless reviewed=yes/no is specified as
    part of the query, all of UniProtKB is searched, which means there
    could be very many hits.

    Args:
        query (str): UniProtKB query string.
        frmt (str; optional): Format for results.
            Can be txt, xml, rdf, fasta. Defaults to txt.
        file_handle (bool): Whether to wrap returned string in StringIO.
            Defaults to False.

    Returns:
        str or None: Data, if any.
    '''
    result = _search(query, frmt=frmt, file_handle=file_handle)
    return result


def number_SP_hits(query):
    '''Search reviewed UniProtKB entries only (Swiss-Prot).

    Accepts standard UniProtKB query syntax.

    Args:
        query (str): UniProtKB query string.

    Returns:
        int: number of hits.
    '''
    result = _search(query, frmt='list', file_handle=False)
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
    result = session.get(UNIPROT_BATCH, params=payload)
    if result.ok:
        if len(result.content) > 0:
            if file:
                return io.StringIO(result.content.decode())
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
    response = session.post(UNIPROT_CONVERT,
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
    response = session.get(UNIPROT_MAP, params=payload)
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


def _search(query, frmt='txt', file_handle=False):
    fmt = frmt.lower()
    _check_format(fmt)
    payload = {'query':query, 'format':fmt}
    result = session.get(UNIPROT_KNOWLEDGEBASE, params=payload)
    if result.ok:
        if len(result.content) > 0:
            if file_handle:
                return io.StringIO(result.text)
            else:
                return str(result.text)
        else:
            raise NoDataError('No data were retrieved.')
    else:
        result.raise_for_status()


def _check_format(fmt):
    allowed_formats = ('html',
                      'tab',
                      'xls',
                      'fasta',
                      'gff',
                      'txt',
                      'xml',
                      'rdf',
                      'list',
                      )
    if fmt not in allowed_formats:
        msg = 'Allowed: {0}\nReceived: {1}'
        raise ValueError(msg.format(allowed_formats, fmt))
