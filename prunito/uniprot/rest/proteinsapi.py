import requests
from ...utils import PROTEINS_API_TAXONOMY


session = requests.Session()


def get_info_on_taxID(taxID):
    """Get details about one taxonomic node as identified by a taxID.

    This includes names, rank, UniProt mnemonics, parent, sibling
    and child nodes, if any. This only accepts a single taxID.
    Use get_info_on_taxIDs() for several.

    Args:
        taxID (str or int): Identifier of the taxonomic node you want info on.
            Can be string or integer.

    Returns:
        dict: Data for taxonomic node.

    Raises:
        NoDataError: If no valid data at all are returned.
    """
    headers = {'Accept': 'application/json'}
    r = session.get('/'.join([PROTEINS_API_TAXONOMY, 'id', str(taxID)]), headers=headers)
    content = r.json()
    if r.status_code == 400 or r.status_code == 404 or r.status_code == 500:
        raise NoDataError(content['errorMessage'])
    else:
        return content


def get_info_on_taxIDs(taxIDs):
    """Get details about one or more taxonomic node(s) as identified by (a) taxID(s).

    This includes names, rank, UniProt mnemonics, parent, sibling
    and child nodes, if any. Currently, any invalid or non-existing taxIDs
    in the query that return errors are silently skipped; we could log this.

    Args:
        taxIDs (Iterable): Iterable of identifiers of the taxonomic nodes you want info on.
            Items in iterable can be string or integer.

    Yields:
        dict: Data for each taxID

    Raises:
        ValueError: If taxIDs are not provided as list or tuple.
        NoDataError: If no valid data at all are returned.
    """
    if isinstance(taxIDs, str):
        raise ValueError('TaxIDs have to be provided as lists or tuples, not strings.')
    else:
        ids_stringified = ','.join([str(item) for item in taxIDs])
        headers = {'Accept': 'application/json'}
        r = session.get('/'.join([PROTEINS_API_TAXONOMY, 'ids', ids_stringified]), headers=headers)
        content = r.json()
        try:
            for node in content['taxonomies']:
                yield node
        except KeyError:
            try:
                raise NoDataError(content['errorMessage'])
            except KeyError:
                raise NoDataError(content['errors'][0]['errorMessage'])



def get_lineage_for_taxID(taxID):
    """Retrieve all ancestor nodes of the given one.

    The data retrieved for each node are taxID and the scientific name.

    Args:
        taxID (str or int): TaxID (NCBI taxonomy identifier)

    Yields:
        dict: Data for each taxID

    Raises:
        NoDataError: If the taxID is invalid or nonexistent.
    """
    headers = {'Accept': 'application/json'}
    r = session.get('/'.join([PROTEINS_API_TAXONOMY, 'lineage', str(taxID)]), headers=headers)
    content = r.json()
    try:
        for node in content['taxonomies']:
            yield node
    except KeyError:
        raise NoDataError(content['errorMessage'][0])


class NoDataError(Exception):
    pass
