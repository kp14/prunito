# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 09:13:33 2015

@author: kp14
"""

import requests
from biocuration.utils import UNIPROT_CONVERT

def convert(path, typ='uniprot', from_='txt', to='xml', encoding='ascii'):
    '''Convert between different data formats using UniProt service.

    Parameters:
        path: path to file containing entry to be converted
        typ: type of the format; default: uniprot
        from: source format
        to: target format
        encding: encoding of the files to be sent; default: ascii

    Returns:
        string
    '''
    payload = {'type': typ,
               'from': from_,
               'to': to
               }
    files = {'data': open(path, 'r', encoding=encoding)}
    response = requests.post(UNIPROT_CONVERT,
                             data=payload,
                             files=files
                             )
    if response.ok:
        return response.text
    else:
        response.raise_for_status()

