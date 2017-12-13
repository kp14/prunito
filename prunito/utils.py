"""Utilities, regexs etc. for the biocuration package.
"""

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

TRANSEQ_BASE = 'http://www.ebi.ac.uk/Tools/services/rest/emboss_transeq'
TRANSEQ_RUN = '/'.join([TRANSEQ_BASE, 'run'])
TRANSEQ_STATUS = '/'.join([TRANSEQ_BASE, 'status'])
TRANSEQ_RESULTS = '/'.join([TRANSEQ_BASE, 'result'])
##############################################################
# General purpose functions
##############################################################

def is_value_in_iterable(val, iterable):
    '''Test whther values is contained in iterable like list or set'''
    return val in iterable


def validate_param(param, val, allowed_vals):
    if not is_value_in_iterable(val, allowed_vals):
        raise ValueError('Wrong parameter {0}!\n'
                         'Allowed:{1}\n'
                         'Got:{2}'.format(param, allowed_vals, val))

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


class NoDataError(Exception):
    pass