# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:25:43 2015

@author: kp14
"""
import os
import pytest

from Bio import SwissProt
from biocuration.uniprot import parse_txt_compatible, parse_txt_atomic


#####################################################################
# Test stuff for parse_txt_compatible
#####################################################################

# Import paths for data files
base_dir = os.path.dirname(__file__)
#filename = 'one_sp_entry.txl'
filename = 'many_sp_entries.txl'
#datafile = os.path.join('SwissProt', filename)
datafile = os.path.join(base_dir, 'SwissProt', filename)

PROPS = ['entry_name',
         'data_class',
         'molecule_type',
         'sequence_length',
         'accessions',
         'created',
         #'sequence_update',
         #'annotation_update',
         'description',
         'gene_name',
         'organism',
         'taxonomy_id',
         'host_organism',
         'host_taxonomy_id',
         'organelle',
         'organism_classification',
         #'references',
         'comments',
         'cross_references',
         'keywords',
         'seqinfo',
         'features',
         'sequence',
         ]

with open(datafile, "r", encoding="ascii") as data:
    my_records = list(parse_txt_compatible(data))


with open(datafile, "r", encoding="ascii") as data:
    biopython_records = list(SwissProt.parse(data))
#####################################################################
# END - Test stuff for parse_txt_compatible
#####################################################################

#####################################################################
# Test stuff for parse_txt_atomic
#####################################################################
EXPECTED_ANNOTATIONS = [
    "P62140: FUNCTION-Dephosphorylates the 'Ser-418' residue of FOXP3 in regulatory T- cells (Treg) from patients with rheumatoid arthritis, thereby inactivating FOXP3 and rendering Treg cells functionally defective - ECO:0000269-PubMed:23396208",
    'P62140: FUNCTION-Protein phosphatase that associates with over 200 regulatory proteins to form highly specific holoenzymes which dephosphorylate hundreds of biological targets - ECO:0000269-None',
    'P62140: FUNCTION-Protein phosphatase (PP1) is essential for cell division, it participates in the regulation of glycogen metabolism, muscle contractility and protein synthesis - ECO:0000269-None',
    'P62140: FUNCTION-Involved in regulation of ionic conductances and long-term synaptic plasticity - ECO:0000269-None',
    'P62140: FUNCTION-Component of the PTW/PP1 phosphatase complex, which plays a role in the control of chromatin structure and cell cycle progression during the transition from mitosis into interphase - ECO:0000269-None',
    'P62140: FUNCTION-In balance with CSNK1D and CSNK1E, determines the circadian period length, through the regulation of the speed and rhythmicity of PER1 and PER2 phosphorylation - ECO:0000269-None',
    'P62140: FUNCTION-May dephosphorylate CSNK1D and CSNK1E - ECO:0000269-None',
    'P62140: CATALYTIC ACTIVITY-[a protein]-serine/threonine phosphate + H(2)O = [a protein]-serine/threonine + phosphate - ECO:0000000-None',
    'P62140: CATALYTIC ACTIVITY-[Myosin light-chain] phosphate + H(2)O = [myosin light-chain] + phosphate - ECO:0000269-PubMed:20516061',
    'P62140: ENZYME REGULATION-Inhibited by the toxins okadaic acid, tautomycin and microcystin Leu-Arg - ECO:0000250-None',
    'P62140: ENZYME REGULATION-The phosphatase activity of the PPP1R15A-PP1 complex toward EIF2S1 is specifically inhibited by Salubrinal, a drug that protects cells from endoplasmic reticulum stress - ECO:0000269-PubMed:15705855',
    'P62140: SUBUNIT-Interacts with RRP1B - ECO:0000269-PubMed:20926688',
    'P62140: SUBUNIT-Part of a complex containing PPP1R15B, PP1 and NCK1/2 - ECO:0000250-UniProtKB:P62141',
    'P62140: SUBUNIT-PP1 comprises a catalytic subunit, PPP1CA, PPP1CB or PPP1CC, which is folded into its native form by inhibitor 2 and glycogen synthetase kinase 3, and then complexed to one or several targeting or regulatory subunits - ECO:0000269-None',
    'P62140: SUBUNIT-The targeting or regulatory subunits determine the substrate specificity of PP1 - ECO:0000269-None',
    'P62140: SUBUNIT-PPP1R12A, PPP1R12B and PPP1R12C mediate binding to myosin - ECO:0000269-None',
    'P62140: SUBUNIT-PPP1R3A (in skeletal muscle), PPP1R3B (in liver), PPP1R3C, PPP1R3D and PPP1R3F (in brain) mediate binding to glycogen - ECO:0000269-None',
    'P62140: SUBUNIT-Component of the MLL5-L complex, at least composed of KMT2E/MLL5, STK38, PPP1CA, PPP1CB, PPP1CC, HCFC1, ACTB and OGT - ECO:0000269-None',
    'P62140: SUBUNIT-Interacts with PPP1R7 and PPP1R12C - ECO:0000269-None',
    'P62140: SUBUNIT-PPP1R15A and PPP1R15B mediate binding to EIF2S1 - ECO:0000269-None',
    'P62140: SUBUNIT-Interacts with PPP1R16B - ECO:0000269-None',
    'P62140: SUBUNIT-Component of the PTW/PP1 phosphatase complex, composed of PPP1R10/PNUTS, TOX4, WDR82, and PPP1CA or PPP1CB or PPP1CC - ECO:0000269-None',
    'P62140: SUBUNIT-Interacts with PPP1R8 - ECO:0000269-None',
    'P62140: SUBUNIT-Interacts with TRIM28; the interaction is weak - ECO:0000269-None',
    'P62140: SUBUNIT-Interacts with PPP1R12A and NUAK1; the interaction is direct - ECO:0000269-None',
    'P62140: SUBUNIT-Interacts with FOXP3 - ECO:0000269-None',
    'P62140: SUBCELLULAR LOCATION-Cytoplasm - ECO:0000269-PubMed:11739654',
    'P62140: SUBCELLULAR LOCATION-Nucleus - ECO:0000269-PubMed:11739654',
    'P62140: SUBCELLULAR LOCATION-Nucleus, nucleoplasm - ECO:0000269-PubMed:11739654',
    'P62140: SUBCELLULAR LOCATION-Nucleus, nucleolus - ECO:0000269-PubMed:11739654',
    'P62140: SUBCELLULAR LOCATION-Nucleus, nucleolus - ECO:0000269-PubMed:20926688',
    'P62140: SUBCELLULAR LOCATION-Highly mobile in cells and can be relocalized through interaction with targeting subunits - ECO:0000000-None',
    'P62140: SUBCELLULAR LOCATION-In the presence of PPP1R8 relocalizes from the nucleus to nuclear speckles - ECO:0000000-None',
    'P62140: INDUCTION-Up-regulated in synovial fluid mononuclear cells and peripheral blood mononuclear cells from patients with rheumatoid arthritis - ECO:0000269-PubMed:23396208',
    'P62140: SIMILARITY-Belongs to the PPP phosphatase family - ECO:0000305-None',
    'P62140: SIMILARITY-PP-1 subfamily - ECO:0000305-None',
    'P62140: INIT_MET-INIT_MET Removed - ECO:0000244-PubMed:22223895',
    'P62140: INIT_MET-INIT_MET Removed - ECO:0000244-PubMed:22814378',
    'P62140: INIT_MET-INIT_MET Removed - ECO:0000269-PubMed:12665801',
    'P62140: CHAIN-CHAIN Serine/threonine-protein phosphatase PP1-beta catalytic subunit. - None',
    'P62140: ACT_SITE-ACT_SITE Proton donor - ECO:0000250-None',
    'P62140: METAL-METAL Manganese 1 - ECO:0000250-None',
    'P62140: METAL-METAL Manganese 1 - ECO:0000250-None',
    'P62140: METAL-METAL Manganese 1 - ECO:0000250-None',
    'P62140: METAL-METAL Manganese 2 - ECO:0000250-None',
    'P62140: METAL-METAL Manganese 2 - ECO:0000250-None',
    'P62140: METAL-METAL Manganese 2 - ECO:0000250-None',
    'P62140: METAL-METAL Manganese 2 - ECO:0000250-None',
    'P62140: MOD_RES-MOD_RES N-acetylalanine - ECO:0000244-PubMed:22223895',
    'P62140: MOD_RES-MOD_RES N-acetylalanine - ECO:0000244-PubMed:22814378',
    'P62140: MOD_RES-MOD_RES N-acetylalanine - ECO:0000269-PubMed:12665801',
    'P62140: MOD_RES-MOD_RES Phosphothreonine - ECO:0000244-PubMed:18669648',
    'P62140: MOD_RES-MOD_RES Phosphothreonine - ECO:0000244-PubMed:23186163',
    'P62140: VARIANT-VARIANT P -> R (probable disease-associated mutation found in patients with rasopathy reminiscent of NSLH) - ECO:0000269-PubMed:27264673',
    'P62140: VARIANT-VARIANT A -> P (probable disease-associated mutation found in patients with rasopathy reminiscent of NSLH) - ECO:0000269-PubMed:27264673',
    'P62140: CONFLICT-CONFLICT L -> P (in Ref. 5; AAV38549) - ECO:0000305-None',
]


#filename = 'one_sp_entry.txl'
filename_atomic = 'atomic.txl'
#datafile = os.path.join('SwissProt', filename)
datafile_atomic = os.path.join(base_dir, 'SwissProt', filename_atomic)

ACTUAL_ANNOTATIONS = []

with open(datafile_atomic, 'r', encoding='ascii') as infile:
    pile = parse_txt_atomic(infile)
    for anno in pile:
        ACTUAL_ANNOTATIONS.append(anno.__str__())
#####################################################################
# END - Test stuff for parse_txt_atomic
#####################################################################


def pytest_generate_tests(metafunc):
    '''Automated generation of tests based on input of multiple UniProt entries.

    This targets the metafunc `test_compare_annotations` in module
    test_uniprotkb_multi.py which just compare two objects for identity, mostly
    strings. Each individual piece of annotation in the provided UniProt entries
    is used to generate its own test which makes pinpointing errors easier and
    saves having to code the test manually. We can also use an arbitray set of
    entries to test.

    ALTERNATIVE PRODUCTS, COFACTOR and BIOPHYSICOCHEMICAL PROPS are right now
    marked as expected failures as my parser
    strips out some white space as opposed to the biopython one.
    '''
    if 'my_annotation' and 'biopy_annotation' in metafunc.fixturenames:
        argvalues = []
        for my_record, record in zip(my_records, biopython_records):
            for prop in PROPS:
                if prop == 'comments':
                    comments = []
                    for my_comment, biopy_comment in zip(getattr(my_record, prop),
                                                         getattr(record, prop)):
                        if my_comment[:8] in ('COFACTOR',
                                              'ALTERNAT',
                                              'BIOPHYSI'):
                            # There is less white space here vs. the biopython
                            # parser, so we expect these comparisons to fail.
                            tple = pytest.mark.xfail((my_comment, biopy_comment))
                            # Add this to the list right away
                            comments.append(tple)
                            # Then go on to construct new tuples with the white
                            # space stripped out which should pass tests.
                            tple = (my_comment.replace(' ', ''),
                                    biopy_comment.replace(' ', ''))
                        else:
                            tple = (my_comment, biopy_comment)
                        comments.append(tple)
                    argvalues.extend(comments)
                # Comment out the following elif block to not check xrefs individually
                elif prop == 'cross_references':
                    xrefs = []
                    for my_xref, biopy_xref in zip(getattr(my_record, prop),
                                                         getattr(record, prop)):
                        tple = (my_xref, biopy_xref)
                        xrefs.append(tple)
                    argvalues.extend(xrefs)
                else:
                    argvalues.append((getattr(my_record, prop),
                                      getattr(record, prop)))
        metafunc.parametrize(['my_annotation', 'biopy_annotation'], argvalues)
    elif 'actual' and 'expected_list' in metafunc.fixturenames:
        argvalues = []
        for anno in ACTUAL_ANNOTATIONS:
            argvalues.append((anno, EXPECTED_ANNOTATIONS))
        metafunc.parametrize(['actual', 'expected_list'], argvalues)
    elif 'expected' and 'actual_list' in metafunc.fixturenames:
        argvalues = []
        for anno in EXPECTED_ANNOTATIONS:
            argvalues.append((anno, ACTUAL_ANNOTATIONS))
        metafunc.parametrize(['expected', 'actual_list'], argvalues)
