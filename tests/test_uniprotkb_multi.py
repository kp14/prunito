#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test for the SwissProt parser on SwissProt files.

Created on Mon May 18 11:59:36 2015

@author: kpichler
"""

import os

from unittest import TestCase

from Bio import SwissProt
import biocuration.uniprotkb.parsing as up

# Some diffs are very long here if tests fail so we set this
# to None to display the whole diff
TestCase.maxDiff = None

# Import paths for data files
base_dir = os.path.dirname(__file__)
#filename = 'one_sp_entry.txl'
filename = 'many_sp_entries.txl'
#datafile = os.path.join('SwissProt', filename)
datafile = os.path.join(base_dir, 'SwissProt', filename)

PROPS = ['entry_name',
         'data_class',
         'sequence_length',
         'accessions',
         'description',
         'gene_name',
         'comments',
         ]

with open(datafile, "r", encoding="ascii") as data:
    my_records = list(up.parse_txt_compatible(data))


with open(datafile, "r", encoding="ascii") as data:
    biopython_records = list(SwissProt.parse(data))


def test_number_of_parsed_entries_is_equal():
    assert len(my_records) == len(biopython_records)

def test_many_entries():
    failures = 0
    for my_record, record in zip(my_records, biopython_records):
        for prop in PROPS:
            my_prop = getattr(my_record, prop)
            biopy_prop = getattr(record, prop)
            try:
                assert  my_prop == biopy_prop
            except AssertionError:
                failures += 1
                print('\n my:{0} <-> biopy: {1}:\n{2}\n{3}'.format(my_record.accessions[0],
                                                               record.accessions[0],
                                                               my_prop,
                                                               biopy_prop))
    assert failures == 0