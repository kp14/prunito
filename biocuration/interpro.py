# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 14:57:52 2015

@author: kp14
"""
import requests
import sys
import time

from .uniprot import search_all
from .utils import InterProXrefs, EBI_HMMER, validate_param
try:
    import venndy.draw as vd
except ImportError:
    sys.exit('Depends on the `venndy` package.')


def search_phmmer(seq=None, seqdb=None, output='json', **kwargs):
    '''Run a protein sequence against a sequence DB using phmmer.

    The input sequence can either be provided in FASTA format or alternatively
    as a retrievable accession number, e.g. a UniprotKB accession.

    Parameters:

    seq: input sequence in FASTA format; string; cannot be used together with `acc`

    acc: accession of sequence to run; cannot be used together with `seq`

    seqdb: sequence database to run search against;
           Possible values: uniprotkb, swissprot, pdb, rp15, rp35, rp55, rp75,
                            uniprotrefprot, pfamseq, qfo

    output: format of output; default: json
            Possible values: json, text, xml, yaml

    hits: humber of results to retrieve; has to be specified as tuple if ints

    return_alignments: set to True of alignments should be retrieved

    e: E value cutoff for returned hits
    ...

    Returns:
    string in chosen format
    '''
    return _search_hmmer(seq=seq, seqdb=seqdb, output=output, **kwargs )


def _search_hmmer(tool='phmmer',
                  seq=None,
                  seqdb=None,
                  output=None,
                  hits=None,
                  return_alignments=False,
                  e=None,
                  domE=None,
                  incE=None,
                  incdomE=None):
    '''Some docs
    '''
    validate_param('seqdb', seqdb, _HMMER_PARAMS['seqdb'])
    validate_param('tool', tool, _HMMER_PARAMS['tool'])
    validate_param('output', output, _HMMER_PARAMS['output'])
    session = requests.Session()
    payload = {'seqdb': seqdb,
               'seq': seq,
              }
    url = '/'.join([EBI_HMMER, tool])
    posted = session.post(url, data=payload, allow_redirects=False)
    if posted.ok:
        time.sleep(3)
        output_values = {'output': output,
                         }
        if hits:
            output_values['range'] = str(hits[0]) + ',' + str(hits[1])
        results = session.get(posted.headers['location'], data=output_values)
        #results = session.get(posted.headers['location'], headers={'Accept': 'application/json'})
        if output == 'json':
            return results.json()
        else:
            return results.content.decode('utf-8')


def draw_signature_overlaps(list_of_signatures, mode='save'):
    '''Represent overlaps in UniprotKB coverage of InterPro xrefs as Venn diagram.

    Parameters:
    list_of_signatures: list of InterPro xref identifiers
    mode: 'save', 'ipython', 'raw'; defaults to save

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

    elif mode == 'save':
        name = '_'.join(list_of_signatures) + '.svg'
        with open(name, 'w') as svg:
            svg.write(rendered_diagram)

    elif mode == 'raw':
        return rendered_diagram


def _get_signature_hit_list(sig):
    '''Retrieve UniProKB hits for a given siganture.

    Parameters:
    sig: InterPro (member) database identifier; string

    Returns:
    set of UniProtKB accessions
    '''
    result = search_all(sig, frmt='list')
    if result:
        result_set = set(result.strip().split('\n'))
    else:
        result_set = set()
        print('Signature returned empty:{}'.format(sig))
    return result_set


_HMMER_PARAMS = {'seqdb': ('pdb',
                           'swissprot',
                           'uniprotkb',
                           'rp15',
                           'rp35',
                           'rp55',
                           'rp75',
                           'uniprotrefprot',
                           'pfamseq',
                           'qfo'
                           ),
                'hmmdb': ('pfam',
                          'gene3d',
                          'tigrfam',
                          'superfam',
                          ),
                'tool': ('phmmer',
                         'hmmsearch',
                         'jackhmmer',
                         'hmmscan',
                         ),
                'output': ('json',
                           'xml',
                           'text',
                           'yaml',
                          ),
                }

if __name__ == '__main__':
    import json
    res = search_phmmer(seq='>Seq\nKLRVLGYHNGEWCEAQTKNGQGWVPSNYITPVNSLENSIDKHSWYHGPVSRNAAEY',
                        seqdb='swissprot', hits=(1,1), output='text')
    print(res)
    #print(json.dumps(res, sort_keys=True, indent=4))
    #draw_signature_overlaps(['PTHR10159', 'PR01908', 'PR01909', 'PR01764', 'PIRSF000939'], mode='save')
