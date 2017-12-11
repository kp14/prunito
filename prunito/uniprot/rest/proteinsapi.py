import requests
from ...utils import PROTEINS_API_TAXONOMY


SESSION = requests.Session()


def tax_node_info(taxIDs, fmt='json'):
    """Get details about one or more taxonomic node(s).

    This includes names, rank, UniProt mnemonics, parent, sibling
    and child nodes, if any.

    Args:
        taxID (iterable): Iterable of identifiers of the taxonomic nodes you want info on.
            Items in iterator can be strings or integers.
        fmt (str): Format information is wanted in. Can be JSON or XML.
            Defaults to JSON.

    Returns:
        Data (JSON dict or XML string)
    """
    fmt = fmt.lower()
    allowed = ('json', 'xml')
    if not fmt in allowed:
        return ValueError('Invalid format: {0}. Choose one of: {1}'.format(fmt,
                                                                           allowed))
    else:
        headers = {'Accept': 'application/{}'.format(fmt)}
        r = SESSION.get('/'.join([PROTEINS_API_TAXONOMY, 'ids', taxIDs]), headers=headers)
        if not r.ok:
            r.raise_for_status()
        else:
            if fmt == 'json':
                return r.json()
            else:
                return r.text
