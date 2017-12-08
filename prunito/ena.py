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

    ProteinIDs (pid) have the format <ID>.<VERSION>, e.g. BAC67592.1
    They are different from ENA accession numbers although those look similar.
    A given ENA entry (as identified by an accession) can contain zero to many
    pid. For retrieval from ENA, versions of pid should be stripped.

    Args:
        pid (str): The protein_ID, e.g. BAC67592.
        fmt (str, optional): Retrieval format. Can be fasta, text, xml. Defaults to fasta.

    Returns:
    CDS nucleotide seq in FastA format (str)
    '''
    try:
        # No checks are currently made to ensure that pid is a pid
        return retrieve(pid, fmt=fmt)
    except ValueError as e:
        print(pid, e)


def retrieve(identifier, fmt='text'):
    '''Retrieve data based on an ENA identifier.

    This is intended for ENA accession numbers and proteinIDs but I guess it
    should work for assembly metadata (e.g. GCA_001521735) or BioProject metadata
    (PRJNA301708), too. Not sure how useful the latter are though.

    Args:
        identifier (str): Identifier to retrieve.
        fmt (str, optional): Retrieval format. Can be fasta, text, xml. Defaults to text.

    Returns:
        Data linked to identifier (str)
    '''
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