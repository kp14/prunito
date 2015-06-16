# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 14:57:52 2015

@author: kp14
"""
import sys
import biocuration.uniprotkb as ukb
from biocuration.utils import InterProXrefs
try:
    import venn.draw as vd
except ImportError:
    sys.exit('Depends on the `venn` package.')


def draw_signature_overlaps(list_of_signatures, mode='save'):
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
    result = ukb.search_all(sig, frmt='list')
    if result:
        result_set = set(result.strip().split('\n'))
    else:
        result_set = set()
        print('Signature returned empty:{}'.format(sig))
    return result_set


if __name__ == '__main__':
    draw_signature_overlaps(['PTHR10159', 'PR01908', 'PR01909', 'PR01764', 'PIRSF000939'], mode='save')