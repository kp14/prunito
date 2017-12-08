# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 10:26:06 2016

@author: kp14
"""

import requests
try:
    from Bio import SeqIO
except ImportError:
    print('No Bio.SeqIO module - functionality relying on parsing EMBL files will be unavailable.')

from . import utils

session = requests.Session()

def retrieve_cds(pid, fmt='fasta'):
    '''Get CDS referenced by a protein_id.

    Parameter:
    pid: pid, string like BAC67592(.1) ecluding the version
    fmt: fasta, text or xml; default to fasta

    Returns:
    CDS nucleotide seq in FastA format
    '''
    try:
        # No checks are currently made to ensure that pid is a pid
        return _retrieve_data(pid, fmt=fmt)
    except ValueError as e:
        print(pid, e)


def _retrieve_data(identifier, fmt='text'):
    '''Retrieva data from ENA via REST.'''
    error_msg = 'display type is either not supported or entry is not found'
    data = identifier + '&display=' + fmt
    result = session.get('/'.join([utils.ENA_DATA, data]))
    if result.ok:
        decoded_result = result.content.decode('utf8')
        if error_msg in decoded_result:
            raise ValueError(error_msg)
        return decoded_result
    else:
        result.raise_for_status()

def embl_cds2fasta(embl_entry, gene_name):
    pass
#from Bio import SeqIO
#
#
#with open('embl.tmp', 'r') as ghandle, open('gluclalpha.fasta', 'w') as out:
#    cds_count = 0
#    for record in SeqIO.parse(ghandle, 'embl'):
#        for ft in record.features:
#            if ft.type == 'CDS':
#                cds_count += 1
#                if 'GluClalpha' in ft.qualifiers['gene']:
#                    out.write('>{0}\n{1}\n'.format(ft.qualifiers['product'], ft.qualifiers['translation'][0]))

def embl2fasta(embl_entry):
    pass
#with open('embl.tmp', 'r') as ghandle, open('embl.fasta', 'w') as out:
#    SeqIO.write(SeqIO.parse(ghandle, 'embl'), out, 'fasta')