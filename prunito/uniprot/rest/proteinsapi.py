import requests
from ...utils import PROTEINS_API_TAXONOMY


SESSION = requests.Session()


def get_info_on_taxID(taxID, fmt='json'):
    """Get details about one or more taxonomic node(s) as identified by a taxID.

    This includes names, rank, UniProt mnemonics, parent, sibling
    and child nodes, if any. This only accepts a single taxID. Use tax_ids_info()
    for several.

    Args:
        taxID (str or int): Identifier of the taxonomic node you want info on.
            Can be string or integer.
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
        r = SESSION.get('/'.join([PROTEINS_API_TAXONOMY, 'id', str(taxID)]), headers=headers)
        if not r.ok:
            r.raise_for_status()
        else:
            if fmt == 'json':
                return r.json()
            else:
                return r.text
