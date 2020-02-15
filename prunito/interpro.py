# -*- coding: utf-8 -*-
"""Module for accessing InterPro via REST-ful web services."""

import requests
import sys
import time

from .uniprot import search
from .utils import EBI_HMMER, NoDataError, WSResponseHmmer

try:
    import venndy as vd
except ImportError:
    import warnings

    msg = (
        "Package venndy not installed. "
        "Drawing InterPro overlaps will not be possible."
    )
    warnings.warn(msg, ImportWarning)

interpro_session = requests.Session()


def search_phmmer(seq, database="swissprot", fmt="json", hits=10, alignments=False):
    """Run a protein sequence against a sequence DB using phmmer.

    The input sequence can either be provided in FASTA format or alternatively
    as a retrievable accession number, e.g. a UniprotKB accession.

    Args:

        seq (str): input sequence in FASTA format.

        database (str): sequence database to run search against;
               Possible values: uniprotrefprot, uniprotkb, swissprot, pdb,
               rp15, rp35, rp55, rp75, ensemblgenomes, ensembl, qfo

        output (str): format of output; default: json
                Possible values: json, text, xml, yaml

        hits (int): humber of results to retrieve; Max: 1000, Default: 10

        alignments (boolean): Whether alignments should be retrieved. Default: False.

    Returns:
        data in selected format (default JSON)
    """
    payload = {
        "seqdb": database,
        "seq": seq,
    }
    url = "/".join([EBI_HMMER, "phmmer"])
    posted = interpro_session.post(url, data=payload, allow_redirects=False)
    if posted.ok:
        try:
            url4results = posted.headers["location"]
        except KeyError:
            raise
        output_values = {
            "output": fmt,
            "range": "1,{}".format(str(hits)),
            "ali": alignments,
        }
        return WSResponseHmmer(interpro_session.get(url4results, params=output_values))
    else:
        posted.raise_for_status()


def draw_signature_overlaps(list_of_signatures, mode="save"):
    """Represent overlaps in UniProtKB coverage of InterPro xrefs as Venn diagram.

    Args:
        list_of_signatures: list of InterPro xref identifiers
        mode: 'save', 'ipython', 'raw'; defaults to save

    Returns:
        Saved rendered SVG Venn diagram; wrapped in SVG() container if mode ipython
        Outputfile name: sig1_sig2_..._sigN.svg
    """
    res_sets = []
    for sig in list_of_signatures:
        res_set = uniprot_hits_for_interpro(sig)
        res_sets.append(res_set)

    rendered_diagram = vd.draw(res_sets, labels=list_of_signatures)

    if mode == "ipython":
        try:
            from IPython.core.display import SVG

            return SVG(data=rendered_diagram)
        except ImportError:
            sys.exit("Drawing depends on IPython.")

    elif mode == "save":
        name = "_".join(list_of_signatures) + ".svg"
        with open(name, "w") as svg:
            svg.write(rendered_diagram)

    elif mode == "raw":
        return rendered_diagram


def uniprot_hits_for_interpro(signature):
    """Retrieve UniProKB hits for a given InterPro signature.

    Args:
        signature (str) : InterPro (member) database identifier

    Returns:
        set: UniProtKB accessions
    """
    try:
        result = search(signature, frmt="list")
    except NoDataError:
        return set()
    else:
        return set(result)
