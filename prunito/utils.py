"""Utilities, regexs etc. for the biocuration package.
"""
import datetime
import io
import re
from enum import Enum

##############################################################
# Base URLs for various resources
##############################################################

UNIPROT_BASE = 'http://www.uniprot.org'

EBI_BASE = 'https://www.ebi.ac.uk'

# Example: http://www.uniprot.org/uniprot/?query=name%3Dtest&sort=score
UNIPROT_KNOWLEDGEBASE = '/'.join([UNIPROT_BASE, 'uniprot'])

# Example: http://www.uniprot.org/keywords/?query=antibiotic
UNIPROT_KEYWORD = '/'.join([UNIPROT_BASE, 'keywords'])

# Example: http://www.uniprot.org/taxonomy/?query=bactus
UNIPROT_TAXONOMY = '/'.join([UNIPROT_BASE, 'taxonomy'])

UNIPROT_BATCH = '/'.join([UNIPROT_BASE, 'batch'])

UNIPROT_CONVERT = '/'.join([UNIPROT_BASE, 'convert'])

UNIPROT_MAP = '/'.join([UNIPROT_BASE, 'mapping'])

UNIPROT_UNIRULE = '/'.join([UNIPROT_BASE, 'unirule'])

PROTEINS_API_BASE = '/'.join([EBI_BASE, 'proteins/api'])

PROTEINS_API_TAXONOMY = '/'.join([PROTEINS_API_BASE, 'taxonomy'])

# Example: http://www.ebi.ac.uk/interpro/search?q=kinase
INTERPRO_SEARCH = '/'.join([EBI_BASE, 'interpro/search'])

# EBI HMMER
EBI_HMMER = '/'.join([EBI_BASE, 'Tools/hmmer/search'])

# Example: http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:0006915
QUICKGO_ID = '/'.join([EBI_BASE, 'QuickGO/GTerm'])

# Example: http://www.ebi.ac.uk/QuickGO/GSearch?q=GO:0006915
QUICKGO_SEARCH = '/'.join([EBI_BASE, 'QuickGO/GSearch'])

# Example: http://www.ebi.ac.uk/ena/data/view/AB107287&display=fasta
ENA_DATA = '/'.join([EBI_BASE, 'ena/data/view'])

EPMC_SEARCH = '/'.join([EBI_BASE, 'europepmc/webservices/rest/search/'])

##############################################################
# General purpose functions
##############################################################

def is_value_in_iterable(val, iterable):
    '''Test whther values is contained in iterable like list or set'''
    return val in iterable

##############################################################
# General purpose regular expressions
##############################################################

# Regex for UniProtKB accession numbers
# This was slightly modified based on test on regex101.com
# Original: [O,P,Q][0-9][A-Z, 0-9]{3}[0-9] | [A-N,R-Z]([0-9][A-Z][A-Z, 0-9]{2}){1,2}[0-9]

UNIPROT_ACCESSION = re.compile(r'[O-Q][0-9][A-Z,0-9]{3}[0-9]|[A-N,R-Z](?:[0-9][A-Z][A-Z,0-9]{2}){1,2}[0-9]')
UNIPROT_ACCESSION_REGEX = UNIPROT_ACCESSION
UNIPROT_ACCESSION_REGEX_STRING = r'[O-Q][0-9][A-Z,0-9]{3}[0-9]|[A-N,R-Z](?:[0-9][A-Z][A-Z,0-9]{2}){1,2}[0-9]'
UNIPROT_EVIDENCE_REGEX = re.compile(r'ECO:[0-9]{7}(?:\|[A-Za-z-]*:[A-Za-z0-9]*)?')
PUBMED_REGEX = re.compile(r'PubMed:[0-9]*')
EC_REGEX = re.compile(r'[1-6]\.[0-9-]*\.[0-9-]*\.[0-9-]*')
EC_INCOMPLETE_REGEX = re.compile(r'[1-6]\.[0-9-]*\.[0-9-]*\.-')
ENA_IDENTIFIER_REGEX = re.compile(r'[A-Z][0-9]{5}|[A-Z]{2}[0-9]{6}|[A-Z]{3}[0-9]{5}')

##############################################################
# General purpose classes
##############################################################

class Borg(object):
    """The Borg design pattern

    A class whose instances share the same state.
    An alternative the Singleton DP, arguably more pythonic.

    """
    _shared = {}

    def __new__(cls,*args,**kwargs):
        self = object.__new__(cls)
        self.__dict__ = cls._shared
        return self


##############################################################
# Enums
##############################################################

class UniProtComments(Enum):
    """Enumeration of comment types used in UniProtKB.

       Order of comments os that used in flat files.
    """
    FUNCTION = 1
    CATALYTIC_ACTIVITY = 2
    COFACTOR = 3
    ENZYME_REGULATION = 4
    BIOPHYSICOCHEMICAL_PROPERTIES = 5
    PATHWAY = 6
    SUBUNIT = 7
    INTERACTION = 8
    SUBCELLULAR_LOCATION = 9
    ALTERNATIVE_PRODUCTS = 10
    TISSUE_SPECIFICITY = 11
    DEVELOPMENTAL_STAGE = 12
    INDUCTION = 13
    DOMAIN = 14
    PTM = 15
    RNA_EDITING = 16
    MASS_SPECTROMETRY = 17
    POLYMORPHISM = 18
    DISEASE = 19
    DISRUPTION_PHENOTYPE = 20
    ALLERGEN = 21
    TOXIC_DOSE = 22
    BIOTECHNOLOGY = 23
    PHARMACEUTICAL = 24
    MISCELLANEOUS = 25
    SIMILARITY = 26
    CAUTION = 27
    SEQUENCE_CAUTION = 28
    WEB_RESOURCE = 29

    def __str__(self):
        return self.name.replace('_', ' ')


class UniProtFeatures(Enum):
    """Enumeration representing feature types used in UniProtKB.

       Order of comments os that used in flat files.
    """
    INIT_MET = 1
    SIGNAL = 2
    PROPEP = 3
    TRANSIT = 4
    CHAIN = 5
    PEPTIDE = 6
    TOPO_DOM = 7
    TRANSMEM = 8
    DOMAIN = 9
    REPEAT = 10
    CA_BIND = 11
    ZN_FING = 12
    DNA_BIND = 13
    NP_BIND = 14
    REGION = 15
    COILED = 16
    MOTIF = 17
    COMPBIAS = 18
    ACT_SITE = 19
    METAL = 20
    BINDING = 21
    SITE = 22
    NON_STD = 23
    MOD_RES = 24
    LIPID = 25
    CARBOHYD = 26
    DISULFID = 27
    CROSSLNK = 28
    VAR_SEQ = 29
    VARIANT = 30
    MUTAGEN = 31
    UNSURE = 32
    CONFLICT = 33
    NON_CONS = 34
    NON_TER = 35
    HELIX = 36
    TURN = 37
    STRAND = 38


class InterProXrefs(Enum):
    '''Enum representing InterPro corssrefs as found in UniProtKB.'''
    Gene3D = '[0-9]\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}'
    HAMAP = 'MF_[0-9]{5}_[A,B]{1}|MF_[0-9]{5}'
    InterPro = 'IP[0-9]{9}'
    PANTHER = 'PTHR[0-9]{5}:SF[0-9]{1,5}|PTHR[0-9]{5}' # 'PTHR[0-9]{5}(?::SF[0-9]{1,5})?'
    PIRSF = 'PIRSF[0-9]{6}'
    PRINTS = 'PR[0-9]{5}'
    PROSITE = 'PS[0-9]{5}'
    Pfam = 'Pf[0-9]{5}'
    ProDom = 'PD[0-9]{6}'
    SMART = 'SM[0-9]{5}'
    SUPFAM = 'SSF[0-9]{5}'
    TIGRFAMs = 'TIGR[0-9]{5}'

    def __init__(self, pattern):
        '''Get precompiled regexes as well:

        hamap_re = InterProXrefs.HAMAP.regex
        match = hamap_re.match('some string')
        '''
        self.regex = re.compile(pattern)

VALID_ID_MAPPINGS = {'ACC',
                    'ID',
                    'UPARC',
                    'NF50',
                    'NF90',
                    'NF100',
                    'GENENAME',
                    'EMBL_ID',
                    'EMBL',
                    'P_ENTREZGENEID',
                    'P_GI',
                    'PIR',
                    'REFSEQ_NT_ID',
                    'P_REFSEQ_AC',
                    'UNIGENE_ID',
                    'PDB_ID',
                    'DISPROT_ID',
                    'BIOGRID_ID',
                    'DIP_ID',
                    'MINT_ID',
                    'STRING_ID',
                    'CHEMBL_ID',
                    'DRUGBANK_ID',
                    'GUIDETOPHARMACOLOGY_ID',
                    'SWISSLIPIDS_ID',
                    'ALLERGOME_ID',
                    'ESTHER_ID',
                    'MEROPS_ID',
                    'MYCOCLAP_ID',
                    'PEROXIBASE_ID',
                    'REBASE_ID',
                    'TCDB_ID',
                    'BIOMUTA_ID',
                    'DMDM_ID',
                    'WORLD_2DPAGE_ID',
                    'DNASU_ID',
                    'ENSEMBL_ID',
                    'ENSEMBL_PRO_ID',
                    'ENSEMBL_TRS_ID',
                    'ENSEMBLGENOME_ID',
                    'ENSEMBLGENOME_PRO_ID',
                    'ENSEMBLGENOME_TRS_ID',
                    'GENEDB_ID',
                    'P_ENTREZGENEID',
                    'KEGG_ID',
                    'PATRIC_ID',
                    'UCSC_ID',
                    'VECTORBASE_ID',
                    'WBPARASITE_ID',
                    'ARACHNOSERVER_ID',
                    'CCDS_ID',
                    'CGD',
                    'CONOSERVER_ID',
                    'DICTYBASE_ID',
                    'ECHOBASE_ID',
                    'ECOGENE_ID',
                    'EUHCVDB_ID',
                    'EUPATHDB_ID',
                    'FLYBASE_ID',
                    'GENECARDS_ID',
                    'GENEREVIEWS_ID',
                    'H_INVDB_ID',
                    'HGNC_ID',
                    'HPA_ID',
                    'LEGIOLIST_ID',
                    'LEPROMA_ID',
                    'MAIZEGDB_ID',
                    'MGI_ID',
                    'MIM_ID',
                    'NEXTPROT_ID',
                    'ORPHANET_ID',
                    'PHARMGKB_ID',
                    'POMBASE_ID',
                    'PSEUDOCAP_ID',
                    'RGD_ID',
                    'SGD_ID',
                    'TAIR_ID',
                    'TUBERCULIST_ID',
                    'WORMBASE_ID',
                    'WORMBASE_PRO_ID',
                    'WORMBASE_TRS_ID',
                    'XENBASE_ID',
                    'ZFIN_ID',
                    'EGGNOG_ID',
                    'GENETREE_ID',
                    'HOGENOM_ID',
                    'HOVERGEN_ID',
                    'KO_ID',
                    'OMA_ID',
                    'ORTHODB_ID',
                    'TREEFAM_ID',
                    'BIOCYC_ID',
                    'REACTOME_ID',
                    'UNIPATHWAY_ID',
                    'CLEANEX_ID',
                    'COLLECTF_ID',
                    'CHITARS_ID',
                    'GENEWIKI_ID',
                    'GENOMERNAI_ID',
                    }


class WSResponse():
    """Wrapper round requests Response objects.

    Allows to add convenience methods in particular for
    different result formats. For example, specifying 'list'
    as the format for a UniProtKB query can then actually return
    a list via a list() method in addition to the text.
    """

    def __init__(self, response):
        self.response = response

    def __getattr__(self, item):
        try:
            return getattr(self.response, item)
        except AttributeError:
            return getattr(self, item)

    def list(self):
        """Return the response text content as a Python list.

        Obviously, this only makes sense if the content is list-like
        in the first place. Like when you download a list of UniProt
        accessions.

        Returns:
            list
        """
        l = self.response.text.strip().split('\n')
        return [x for x in l if x]

    def as_file_object(self):
        """Return string content wrapped in ioStringIO.
        Returns:
            StringIO instance
        """
        return io.StringIO(self.text)

    def __repr__(self):
        return '<{0}(requests.get({1}))'.format(self.__class__.__qualname__,
                                                self.response.url)

    def __str__(self):
        return '{0} instance based on request\n{1}'.format(self.__class__.__qualname__,
                                                           self.response.url)


class WSResponseUniprot(WSResponse):
    """Results from uniprot.org REST API.

    In addition to what WSResponse provides, a few custom methods are
    implemented for getting the UniProtKB release date, the release
    number or the number of hits in a query. Although search results
    are basically sequence-like (and __len__ is implemented) making
    all the sequence methods available would necessitate parsing the
    results first which, given the various formats, is why it is not
    implemented for now.

    As WSResponse wraps requests.Response, attributes not found in the
    former will be passed on to the latter. So, the results can be
    accessed as text via the attribute `text`, the URL as `url` etc.
    """

    def __init__(self, response):
        super().__init__(response)

    def release(self):
        """Return UniProt release, eg. 2017_09."""
        return self.response.headers['x-uniprot-release']

    def date(self, as_text=True):
        """Return date of UniProt release.

        Args:
            as_text (bool): Whether to return the date as a string.
                If false, convert to datetime object.

        Returns:
            str or datetime object
        """
        date_in_header = self.response.headers['last-modified']
        if as_text:
            return date_in_header
        else:
            return _convert_date_string(date_in_header)

    def size(self):
        """Number of query hits."""
        return self.__len__()

    def __len__(self):
        return int(self.response.headers['x-total-results'])


class WSResponseEPMC(WSResponse):
    """Results from europePMC REST API.

    For europePMC, results are always retrieved as JSON. If
    there are any hits at all, these are provided as a list
    in the JSON which is why we can easily iterate over or
    slice them. Note that other sequence operations might
    not be supported.
    """

    def __init__(self, response):
        super().__init__(response)

    def size(self):
        """Number of query hits."""
        return self.__len__()

    def __len__(self):
        return self.response.json()['hitCount']

    def __iter__(self):
        return self.response.json()['resultList']['result'].__iter__()

    def __getitem__(self, item):
        return self.response.json()['resultList']['result'].__getitem__(item)


class WSResponseTax(WSResponse):
    """Taxonomy nodes via Proteins API.

    Query results are always retrieved as JSON. If
    there are any hits at all, these are provided as a list
    in the JSON which is why we can easily iterate over or
    slice them. Note that other sequence operations might
    not be supported.
    """

    def __init__(self, response):
        super().__init__(response)

    def size(self):
        """Number of tax nodes contained in results."""
        return self.__len__()

    def __len__(self):
        return len(self.response.json()['taxonomies'])

    def __iter__(self):
        return self.response.json()['taxonomies'].__iter__()

    def __getitem__(self, item):
        return self.response.json()['taxonomies'].__getitem__(item)


class NoDataError(Exception):

    def __init__(self, status_code):
        self.status_code = status_code

    def __str__(self):
        msg = ('No data were retrieved. '
               'Query returned code: {}. '
               'If the code was 200 the query probably had no hits.')
        return msg.format(str(self.status_code))


class ExcessiveDataError(Exception):

    def __init__(self, limit, actual):
        self.limit = limit
        self.actual = actual

    def __str__(self):
        msg = ('Number of hits exceeds limit. '
               'Limit: {0}. Actual hits: {1}. '
               'Please adjust limit or query.')
        return msg.format(str(self.limit), str(self.actual))


def _convert_date_string(date_string):
    """Try and convert date into datetime object"""
    full_date = '%a, %d %b %Y %H:%M:%S GMT'
    simple_date = '%d %b %Y'
    try:
        return datetime.datetime.strptime(date_string, full_date)
    except ValueError:
        return datetime.datetime.strptime(date_string[5:16], simple_date)