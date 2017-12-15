# -*- coding: utf-8 -*-
"""Module for accessing ENA data via REST-ful web services."""
import re
import time
import requests


from .utils import TRANSEQ_RUN, TRANSEQ_STATUS, TRANSEQ_RESULTS, ENA_DATA, NoDataError, ENA_IDENTIFIER_REGEX

session = requests.Session()


def retrieve(identifier, fmt='fasta'):
    '''Retrieve data based on an ENA identifier.

    This is intended for ENA accession numbers and proteinIDs but I guess it
    should work for assembly metadata (e.g. GCA_001521735) or BioProject metadata
    (PRJNA301708), too. Not sure how useful the latter are though.

    Args:
        identifier (str): Identifier to retrieve.
        fmt (str, optional): Retrieval format. Can be fasta, text, xml. Defaults to fasta.

    Returns:
        Data linked to identifier (str)

    Raises:
        NoDataError: If no data are returned
    '''
    match = re.match(ENA_IDENTIFIER_REGEX, identifier)
    if not match or len(match.group()) is not len(identifier):
        raise ValueError('Identifier seems to have wrong format. Provide only 1 identifier.')
    error_msg = 'display type is either not supported or entry is not found'
    data = identifier + '&display=' + fmt
    result = session.get('/'.join([ENA_DATA, data]))
    if result.ok:
        decoded_result = result.content.decode('utf8')
        if error_msg in decoded_result:
            raise NoDataError(error_msg)
        return decoded_result
    else:
        result.raise_for_status()


def translate(seq, frame='1', trim=False):
    """Translate a nucleotide sequence.

    This uses the EMBOSS transeq tool, so in principle its documentation
    applies here, too. Not all possible parameters are currently exposed.

    Args:
        seq (str): Sequence to be translated.
            Transeq allows various formats among which FASTA and EMBL
            are most relevant here.
        frame (str): which frames to translate. Defaults to 6, meaning all.
            Other options are: 1, 2, 3, -1, -2, -3, F (all forward) and R.
        trim (bool): Whether transeq should remove trailing * or X.

    Returns:
         str: Protein sequences in FASTA format.

    Raises:
        NoDataError: If translation fails or any errors using webservice.
    """
    seconds_before_checking_again = 5
    job_id = None
    status = None
    data = {'email': 'kp14git@hotmail.com',
            'trim': trim,
            'frame': frame,
            'sequence': seq}
    job = session.post(TRANSEQ_RUN, data=data)
    if job.ok:
        job_id = job.text
    else:
        job.raise_for_status()
    while not status == 'FINISHED':
        print(status)
        time.sleep(seconds_before_checking_again)
        status_url = '/'.join([TRANSEQ_STATUS, job_id])
        print(status_url)
        r = session.get(status_url)
        status = r.text
        if status in ('ERROR', 'FAILURE'):
            raise NoDataError('Transeq returned: {}'.format(status))
    result_url = '/'.join([TRANSEQ_RESULTS, job_id, 'out'])
    result = session.get(result_url)
    return result.text
