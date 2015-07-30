"""Utilities, regexs etc. for the biocuration package.
"""

import re
from enum import Enum

##############################################################
# Base URLs for various resources
##############################################################

UNIPROT_BASE = 'http://www.uniprot.org'

# Example: http://www.uniprot.org/uniprot/?query=name%3Dtest&sort=score
UNIPROT_KNOWLEDGEBASE = '/'.join([UNIPROT_BASE, 'uniprot'])

# Example: http://www.uniprot.org/keywords/?query=antibiotic
UNIPROT_KEYWORD = '/'.join([UNIPROT_BASE, 'keywords'])

# Example: http://www.uniprot.org/taxonomy/?query=bactus
UNIPROT_TAXONOMY = '/'.join([UNIPROT_BASE, 'taxonomy'])

UNIPROT_BATCH = '/'.join([UNIPROT_BASE, 'batch'])

UNIPROT_CONVERT = '/'.join([UNIPROT_BASE, 'convert'])

# Example: http://www.ebi.ac.uk/interpro/search?q=kinase
INTERPRO_SEARCH = "http://www.ebi.ac.uk/interpro/search"

# EBI HMMER
EBI_HMMER = 'http://www.ebi.ac.uk/Tools/hmmer/search/'

# Example: http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:0006915
QUICKGO_ID = "http://http://www.ebi.ac.uk/QuickGO/GTerm"

# Example: http://www.ebi.ac.uk/QuickGO/GSearch?q=GO:0006915
QUICKGO_SEARCH = "http://http://www.ebi.ac.uk/QuickGO/GSearch"


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
UNIPROT_ACCESSION = re.compile("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")


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
        self.regex = re.compile(pattern)

