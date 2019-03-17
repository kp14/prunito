import io
import itertools
import warnings
from collections import defaultdict
import requests
from ..parsers.parser_knowledgebase_txt import parse_txt
from ...utils import (UNIPROT_KNOWLEDGEBASE,
                      UNIPROT_BATCH,
                      UNIPROT_CONVERT,
                      UNIPROT_MAP,
                      GFF_COLUMNS,
                      VALID_ID_MAPPINGS,
                      NoDataError,
                      WSResponse,
                      _convert_date_string,
                      )
try:
    import pandas as pd
except ImportError:
    warnings.warn('Pandas not available. Exporting to dataframes will not be possible.')


session = requests.Session()


class WSResponseUniprot(WSResponse):
    """Results from uniprot.org REST API.

    In addition to what WSResponse provides, a few custom methods are
    implemented for getting the UniProtKB release date, the release
    number or the number of hits in a query. As search results
    are basically sequence-like--we expect one or more hits after all--
    __len__ is implemented as well as __iter__ for formats we can easily
    parse. For tabular data, export as a dataframe is possible if pandas
    is available.

    As WSResponse wraps requests.Response, attributes not found in the
    former will be passed on to the latter. So, the results can be
    accessed as text via the attribute `text`, the URL as `url` etc.
    """

    def __init__(self, response):
        super().__init__(response)
        self._iter_type = None  # populated by search()

    def release(self):
        """Return UniProt release, eg. 2017_09."""
        return self.response.headers['x-uniprot-release']

    def date(self, as_datetime=False):
        """Return date of UniProt release.

        Args:
            as_datetime (bool): Whether to return the date as a datetime object.
                If true, convert to datetime object.

        Returns:
            str or datetime object
        """
        date_in_header = self.response.headers['last-modified']
        if not as_datetime:
            return date_in_header
        else:
            return _convert_date_string(date_in_header)

    def size(self):
        """Number of query hits."""
        return self.__len__()

    def as_dataframe(self):
        """Return tabular results as dataframe.

        This is only relevant for formats list, tab and gff.

        Returns:
            pandas.Dataframe
        """
        if self._iter_type == 'list':
            return pd.read_csv(self.as_file_object(), sep=',', header=None)
        elif self._iter_type == 'tab':
            return pd.read_csv(self.as_file_object(), sep='\t', header=0)
        elif self._iter_type == 'gff':
            temp = io.StringIO()
            for line in self.as_file_object():
                if not line.startswith('#'):
                    temp.write(line)
            temp.seek(0)
            return pd.read_csv(temp, sep='\t', header=None, names=GFF_COLUMNS)
        else:
            msg = 'Dataframes not supported for non-tabular data: {}'.format(self._iter_type)
            raise NotImplementedError(msg)

    def __len__(self):
        return int(self.response.headers['x-total-results'])

    def __iter__(self):
        if self._iter_type == 'txt':
            for entry in parse_txt(self):
                yield entry
        elif self._iter_type in ('tab', 'list', 'gff'):
            for line in self.as_file_object():
                yield line.strip()
        else:
            msg = 'Iteration not implemented for format: '.format(self._iter_type)
            raise NotImplementedError(msg)


class WSResponseUniprotMapping(WSResponseUniprot):
    """Results from calling UniProt identifier mapping."""

    def __init__(self, response):
        super().__init__(response)

    def as_dict(self):
        """Return mapping results as a dictionary."""
        mapped = defaultdict(list)
        lines = list(self)
        for line in lines[1:]:  # skip header
            try:
                source, target = line.strip().split('\t')
            except ValueError:
                raise
            else:
                mapped[source].append(target)
        return mapped

    def target_ids(self):
        """Return mapped-to target IDs as a list.

        Elements in the list will be unique.

        Returns:
            list (str)
        """
        return list(set(itertools.chain.from_iterable(self.as_dict().values())))


def current_release():
    """Get current public release of UniProtKB.

    Each response returned by a UniProt query also contains this info.
    Releases are specified as <year>_<release>, e.g. 2016_02.
    There are 11 or 12 releases per calendar year.

    Returns:
        str
    """
    # Release number contained in request header
    # So we retrieve a random human entry and look up the value
    r = search_reviewed('*', random='yes', frmt='list')
    return r.release()


def current_release_date(as_datetime=False):
    """Get date of current public release.

    The current release date is available from every UniProt
    query result. So this is just a convenience function.

    Args:
        as_datetime (bool): Whether to return the date as a datetime object.
            If true, convert to datetime object.

    Returns:
        str or datetime object
    """
    r = search_reviewed('*', random='yes', frmt='list')
    return r.date(as_text=as_datetime)


def search_reviewed(query, **kwargs):
    '''Search reviewed UniProtKB only.

    See documentation for search() for details.

    Args:
        query (str): UniProtKB query string.
        **kwargs: Various parameters for the REST API.

    Returns:
        WSResponseUniprot object
    '''
    if 'reviewed:yes' not in query.lower():
        query += ' AND reviewed:yes'
    return search(query, **kwargs)


def search_unreviewed(query, **kwargs):
    '''Search unreviewed UniProtKB only.

    See documentation for search() for details.

    Args:
        query (str): UniProtKB query string.
        **kwargs: Various parameters for the REST API.

    Returns:
        WSResponseUniprot object
    '''
    if 'reviewed:no' not in query.lower():
        query += ' AND reviewed:no'
    return search(query, **kwargs)


def search(query, frmt='txt', limit=2000, **kwargs):
    '''Search UniProtKB.

        Accepts standard UniProtKB query syntax, so queries can be tested
        on the uniprot.org website. Unless reviewed=yes/no is specified as
        part of the query, all of UniProtKB is searched, which means there
        could be very many hits.

        Additional parameters which can be used:

        * random: yes/no. Retrieve one random entry from query set.
        * columns: comma-separated list of UniProtKB field names for use with
          tab-separated format. There must not be any empty spaces around the
          commas. Example values: citation, clusters, comments,
          ec, comment(FUNCTION), features, feature(ACTIVE SITE). For a
          complete overview, go here: https://www.uniprot.org/help/uniprotkb_column_names
        * include: yes/no. Include isoform sequences when format=fasta.
        * compress: yes/no. Return results gzipped.
        * offset: integer. Offset of the first result, typically used together with
          the limit parameter.

        Args:
            query (str): UniProtKB query string.
            frmt (str; optional): Format for results.
                Can be txt, xml, rdf, fasta, tab, list. Defaults to txt.
            limit (int): Maximum number of hits to retrieve. Default is 2000.
            **kwargs: The REST API accepts a number of further parameters. See
                doc string above for details.

        Returns:
            WSResponseUniprot object

        Raises:
            NoDataError: If no results are returned.
        '''
    fmt = frmt.lower()
    _check_format(fmt)
    payload = {'query':query, 'format':fmt, 'limit':limit}
    payload.update(**kwargs)
    with session.get(UNIPROT_KNOWLEDGEBASE, params=payload, stream=True) as r:
        result = WSResponseUniprot(r)
        result._iter_type = fmt
        if result.ok:
            if result.size() > limit:
                msg = ('Partial dataset retrieved. Size: {0}. Retrieved: {1}.\n'
                      'Consider increasing the limit and/or using offset.')
                print(msg.format(result.size(), limit))
            _ = result.content
            return result
        else:
            raise NoDataError(result.status_code)


def retrieve_batch(ac_list, frmt='txt'):
    '''Batch retrieval of uniProtKB entries.

    Args:
        ac_list (list): UniProtKB accessions
        frmt (str; optional): Format for results.
            Can be txt, xml, rdf, fasta. Defaults to txt.

    Returns:
        WSResponseUniprot

    Raises:
        NoDataError
    '''
    id_list = ' '.join(ac_list)
    payload = {'from': 'ACC',
               'to': 'ACC',
               'format': frmt,
               'query': id_list,
               }
    result = WSResponseUniprot(session.get(UNIPROT_BATCH, params=payload))
    result._iter_type = frmt
    if result.ok:
        if len(result.content) > 0:
            return result
        else:
            raise NoDataError(result.status_code)
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
        WSResponseUniprot
    '''
    payload = {'type': 'uniprot',
               'from': source_fmt,
               'to': target_fmt,
               }
    files = {'data': open(path, 'r', encoding=encoding)}
    response = WSResponseUniprot(session.post(UNIPROT_CONVERT,
                                              data=payload,
                                              files=files
                                             ))
    if response.ok:
        return response
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

    See https://www.uniprot.org/help/api_idmapping for details.
    Note: The response.url field contains the URL from which to download
    the mapping, e.g. http://www.uniprot.org/mapping/M20160504763V34ZKX0.tab

    Here is a list of some mapping source/targets (may be out of ate):
    data set - value - direction
    UniProtKB AC/ID	ACC+ID	from
    UniProtKB AC	ACC	both
    UniProtKB ID	ID	both
    UniParc	UPARC	both
    UniRef50	NF50	both
    UniRef90	NF90	both
    UniRef100	NF100	both
    Gene name	GENENAME	both
    CRC64	CRC64	both
    EMBL/GenBank/DDBJ	EMBL_ID	both
    EMBL/GenBank/DDBJ CDS	EMBL	both
    Entrez Gene (GeneID)	P_ENTREZGENEID	both
    RefSeq Nucleotide	REFSEQ_NT_ID	both
    RefSeq Protein	P_REFSEQ_AC	both
    PDB	PDB_ID	both
    ChEMBL	CHEMBL_ID	both
    Allergome	ALLERGOME_ID	both
    MEROPS	MEROPS_ID	both
    Ensembl Protein	ENSEMBL_PRO_ID	both
    Ensembl Transcript	ENSEMBL_TRS_ID	both
    Ensembl Genomes	ENSEMBLGENOME_ID	both
    Ensembl Genomes Protein	ENSEMBLGENOME_PRO_ID	both
    Ensembl Genomes Transcript	ENSEMBLGENOME_TRS_ID	both
    GeneID (Entrez Gene)	P_ENTREZGENEID	both
    KEGG	KEGG_ID	both
    UCSC	UCSC_ID	both
    VectorBase	VECTORBASE_ID	both
    Araport	ARAPORT_ID	both
    CCDS	CCDS_ID	both
    dictyBase	DICTYBASE_ID	both
    FlyBase	FLYBASE_ID	both
    HGNC	HGNC_ID	both
    HPA	HPA_ID	both
    MGI	MGI_ID	both
    MIM	MIM_ID	both
    RGD	RGD_ID	both
    SGD	SGD_ID	both
    WormBase	WORMBASE_ID	both
    WormBase Protein	WORMBASE_PRO_ID	both
    WormBase Transcript	WORMBASE_TRS_ID	both
    Xenbase	XENBASE_ID	both
    ZFIN	ZFIN_ID	both
    GeneTree	GENETREE_ID	both
    OMA	OMA_ID	both
    OrthoDB	ORTHODB_ID	both
    TreeFam	TREEFAM_ID	both
    BioCyc	BIOCYC_ID	both
    Reactome	REACTOME_ID	both

    Args:
        id_list (list): Identifiers to be mapped. Identifiers should be strings.
        source_fmt (str): Format of the provided identifiers. See UniProt help
            for allowed values.
        target_fmt (str): Desired format of the identifiers. See UniProt help
            for allowed values.

    Returns:
        WSResponseUniprot: With map attribute, a dict with source IDs as keys
            and a list of target IDs a values.

    Raises:
        ValueError: If invalid identifier formats are used.
        ValueError: If not at least source or target is UniProt accession/ID
    '''
    _validate_mapping_partners(source_fmt, target_fmt)
    id_list = ' '.join(id_list)
    payload = {'from': source_fmt.upper(),
               'to': target_fmt.upper(),
               'format': 'tab',
               'query': id_list,
               }
    response = WSResponseUniprotMapping(session.get(UNIPROT_MAP, params=payload))
    response._iter_type = 'tab'
    if response.ok:
        return response
    else:
        response.raise_for_status()


def _validate_mapping_partners(source, target):
    if source.upper() not in VALID_ID_MAPPINGS or target.upper() not in VALID_ID_MAPPINGS:
        msg = 'Invalid mapping source(s): {0}, {1}. Use one of:\n{2}'
        raise ValueError(msg.format(source.upper(), target.upper(), ', '.join(VALID_ID_MAPPINGS)))
    if source.upper() not in ('ACC', 'ID') and target.upper() not in ('ACC', 'ID'):
        raise ValueError('One of source or target format has to be UniProt ACC or ID.')


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
