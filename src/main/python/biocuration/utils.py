"""Utilities, regexs etc. for the biocuration package.
"""

import re
from enum import Enum

##############################################################
# Base URLs for various resources
##############################################################

# Example: http://www.uniprot.org/uniprot/?query=name%3Dtest&sort=score
UNIPROT_KNOWLEDGEBASE = "http://www.uniprot.org/uniprot/"

# Example: http://www.uniprot.org/keywords/?query=antibiotic
UNIPROT_KEYWORD = "http://www.uniprot.org/keywords/"

# Example: http://www.uniprot.org/taxonomy/?query=bactus
UNIPROT_TAXONOMY = "http://www.uniprot.org/taxonomy/"

# Example: http://www.ebi.ac.uk/interpro/search?q=kinase
INTERPRO_SEARCH = "http://www.ebi.ac.uk/interpro/search"

# Example: http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:0006915
QUICKGO_ID = "http://http://www.ebi.ac.uk/QuickGO/GTerm"

# Example: http://www.ebi.ac.uk/QuickGO/GSearch?q=GO:0006915
QUICKGO_SEARCH = "http://http://www.ebi.ac.uk/QuickGO/GSearch"


##############################################################
# General purpose regular expressions
##############################################################

# Regex for UniProtKB accession numbers
UNIPROT_ACCESSION = re.compile("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")

# Regex for InterPro entry IDs
INTERPRO_ID = re.compile("IP[0-9]{9}")

# Regex for Pfam signature IDs
PFAM_SIG = re.compile("Pf[0-9]{5}")

# Regex for Hamap signature IDs
# does not work with appended _A/B yet
HAMAP_SIG = re.compile("MF_[0-9]{5}")


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
    CATALYTIC ACTIVITY = 2
    COFACTOR = 3
    ENZYME REGULATION = 4
    BIOPHYSICOCHEMICAL PROPERTIES = 5
    PATHWAY = 6
    SUBUNIT = 7
    INTERACTION = 8
    SUBCELLULAR LOCATION = 9
    ALTERNATIVE PRODUCTS = 10
    TISSUE SPECIFICITY = 11
    DEVELOPMENTAL STAGE = 12
    INDUCTION = 13
    DOMAIN = 14
    PTM = 15
    RNA EDITING = 16
    MASS SPECTROMETRY = 17
    POLYMORPHISM = 18
    DISEASE = 19
    DISRUPTION PHENOTYPE = 20
    ALLERGEN = 21
    TOXIC DOSE = 22
    BIOTECHNOLOGY = 23
    PHARMACEUTICAL = 24
    MISCELLANEOUS = 25
    SIMILARITY = 26
    CAUTION = 27
    SEQUENCE CAUTION = 28
    WEB RESOURCE = 29


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
