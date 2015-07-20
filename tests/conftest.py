# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:25:43 2015

@author: kp14
"""
import os

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

def pytest_generate_tests(metafunc):
    if 'my_annotation' and 'biopy_annotation' in metafunc.fixturenames:
        argvalues = []
        for my_record, record in zip(my_records, biopython_records):
            for prop in PROPS:
                argvalues.append((getattr(my_record, prop), getattr(record, prop)))
        metafunc.parametrize(['my_annotation', 'biopy_annotation'], argvalues)
