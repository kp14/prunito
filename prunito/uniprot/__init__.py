# -*- coding: utf-8 -*-
"""Package for accessing UniProt data via REST-ful web services."""

from .parser_knowledgebase_txt import parse_txt, parse_txt_compatible, parse_txt_atomic
from .parser_unirule import parse_rules
from .rest_api import (current_release,
                       search_reviewed,
                       search_unreviewed,
                       search_all,
                       number_SP_hits,
                       retrieve_batch,
                       convert,
                       map_id,
                       )
from .rest.proteinsapi import (get_info_on_taxID,
                                get_info_on_taxIDs
                               )
from .scraping import get_lineage