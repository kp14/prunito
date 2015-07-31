# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 14:57:52 2015

@author: kp14
"""
import requests
import sys
import time

import biocuration.uniprotkb.searching as us
from biocuration.utils import InterProXrefs, EBI_HMMER
try:
    import venn.draw as vd
except ImportError:
    sys.exit('Depends on the `venn` package.')


def search_phmmer(seq=None, seqdb=None, fmt='json'):
    return _search_hmmer(seq=seq, seqdb=seqdb, fmt=fmt )


def _search_hmmer(seq=None, seqdb=None, tool='phmmer', fmt='json'):
    session = requests.Session()
    payload = {'seqdb': seqdb,
               'seq': seq,
              }
    url = '/'.join([EBI_HMMER, tool])
    posted = session.post(url, data=payload, allow_redirects=False)
    if posted.ok:
        time.sleep(5)
        output_values = {'output': 'json',
                         'range': '1,10', # return only 10 hits
                         }
        results = session.get(posted.headers['location'], data=output_values)
        #results = session.get(posted.headers['location'], headers={'Accept': 'application/json'})
        return results.json()


def draw_signature_overlaps(list_of_signatures, mode='save'):
    '''Represent overlaps in UniprotKB coverage of InterPro xrefs as Venn diagram.

    Parameters:
    list_of_signatures: list of InterPro xref identifiers
    mode: 'save' or 'ipython'; defaults to save

    Returns:
    Saved rendered SVG Venn diagram; wrapped in SVG() container if mode ipython
    Outputfile name: sig1_sig2_..._sigN.svg
    '''
    res_sets = []
    for sig in list_of_signatures:
        res_set = _get_signature_hit_list(sig)
        res_sets.append(res_set)

    rendered_diagram = vd.draw(res_sets, labels=list_of_signatures)

    if mode == 'ipython':
        try:
            from IPython.display import SVG
            return SVG(data=rendered_diagram)
        except ImportError:
            sys.exit('Drawing depends on IPython.')

    if mode == 'save':
        name = '_'.join(list_of_signatures) + '.svg'
        with open(name, 'w') as svg:
            svg.write(rendered_diagram)


def _get_signature_hit_list(sig):
    '''Retrieve UniProKB hits for a given siganture.

    Parameters:
    sig: InterPro (member) database identifier; string

    Returns:
    set of UniProtKB accessions
    '''
    result = us.search_all(sig, frmt='list')
    if result:
        result_set = set(result.strip().split('\n'))
    else:
        result_set = set()
        print('Signature returned empty:{}'.format(sig))
    return result_set


if __name__ == '__main__':
    import json
    res = search_phmmer(seq='>Seq\nKLRVLGYHNGEWCEAQTKNGQGWVPSNYITPVNSLENSIDKHSWYHGPVSRNAAEY',
                        seqdb='swissprot')
    print(json.dumps(res, sort_keys=True, indent=4))
    #draw_signature_overlaps(['PTHR10159', 'PR01908', 'PR01909', 'PR01764', 'PIRSF000939'], mode='save')