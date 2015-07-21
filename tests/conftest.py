# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:25:43 2015

@author: kp14
"""
import os
import pytest

from Bio import SwissProt
import biocuration.uniprotkb.parsing as up

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
    my_records = list(up.parse_txt_compatible(data))


with open(datafile, "r", encoding="ascii") as data:
    biopython_records = list(SwissProt.parse(data))

def pytest_generate_tests(metafunc):
    '''Automated generation of tests based on input of multiple UniProt entries.

    This targets the metafunc `test_compare_annotations` in module
    test_uniprotkb_multi.py which just compare two objects for identity, mostly
    strings. Each individual piece of annotation in the provided UniProt entries
    is used to generate its own test which makes pinpointing errors easier and
    saves having to code the test manually. We can also use an arbitray set of
    entries to test.

    ALTERNATIVE PRODUCTS are right now marked as expected failures as my parser
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
