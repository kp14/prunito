"""Access (some) EBI webservices.

This module provides concrete implementations for accessing EBI
webservices. It also provides a base class that can be subclassed.
"""
import time
import warnings
from lxml import etree
import requests
from .utils import NoDataError, WSResponse

session = requests.Session()

def set_email(email):
    """Provide an email address for EBI webservices."""
    EBIWebService.email = email


class EBIWebService():
    """Base class for generic EBI web service access."""

    email = ''

    def __init__(self, service):
        self._service = service
        self._base_url = 'http://www.ebi.ac.uk/Tools/services/rest/'
        self.url = self._base_url + service
        self._run_url = '/'.join([self.url, 'run'])
        self._status_url = '/'.join([self.url, 'status'])
        self._result_url = '/'.join([self.url, 'result'])
        self._seconds_before_checking_again = 3
        self._finished = False
        self.job_id = None

    def get_parameters(self):
        """Retrieve a list of all parameters for the web service.

        To learn more about possible values for any given parameter,
        use the get_parameter_details() method.

        Returns:
            list of str: Parameters.
        """
        tree = etree.parse('/'.join([self.url, 'parameters']))
        params = []
        for ele in tree.iterfind('./id'):
            params.append(ele.text)
        return params

    def get_parameter_details(self, param):
        """Get details for parameter.

        Details include the name, description, type and possible values.

        Args:
            param (str): A parameter for the web service.

        Returns:
            dict: Parameter details.
        """
        tree = etree.parse('/'.join([self.url, 'parameterdetails', param]))
        root = tree.getroot()
        details = {'values': []}
        for ele in root.iterchildren():
            if ele.tag in {'name', 'description', 'type'}:
                details[ele.tag] = ele.text.replace('\n\t', '')
        for ele in tree.findall('.//value'):
            l = []
            for e in ele.iterchildren('value', 'label'):
                l.append(e.text)
            if l:
                details['values'].append(tuple(l))
        return details

    def _reset(self):
        """Reset parameters before next call."""
        self._finished = False
        self.job_id = None

    def _start_job(self, data):
        """Start web service job.

        Args:
            data (dict): Parameters for web service run.

        Raises:
            HTTPError: Passed on from request.
        """
        self._reset()
        data_all = {'email': EBIWebService.email,}
        data_all.update(data)
        job = session.post(self._run_url, data=data_all)
        if job.ok:
            self.job_id = job.text
        else:
            job.raise_for_status()

    def _is_job_finished(self):
        """Poll status of web service job.

        Returns:
            bool

        Raises:
            NoDataError: If job fails, returns an error or is gone.
        """
        status_url = '/'.join([self._status_url, self.job_id])
        r = session.get(status_url)
        status = r.text
        if status == 'FINISHED':
            self._finished =  True
        elif status == 'RUNNING':
            self._finished = False
        elif status in ('ERROR', 'FAILURE', 'NOT_FOUND'):
            raise NoDataError('{0} returned: {1}'.format(self._service, status))

    def _retrieve_job_results(self):
        """Return protein sequence.

        Returns:
            WSResponse
        """
        result_url = '/'.join([self._result_url, self.job_id, 'out'])
        result = session.get(result_url)
        return WSResponse(result)

    def _run(self, data):
        if not EBIWebService.email:
            print('No email address for EBI web service set.\n'
                  'This is required for the services to work.\n'
                  'Use set_email() to provide an email address.')
        else:
            self._start_job(data)
            while not self._finished:
                time.sleep(self._seconds_before_checking_again)
                self._is_job_finished()
            return self._retrieve_job_results()

    def __call__(self, *args, **kwargs):
        """Override to provide data and signature."""
        data = {}
        return self._run(data)


class EmbossTranseq(EBIWebService):
    """Translate a nucleotide sequence to protein.

    This uses the EMBOSS transeq tool, so in principle its documentation
    applies here, too. Not all possible parameters are currently exposed.

    Documentation for the EMBOSS transeq tool can be found here:
    http://emboss.sourceforge.net/apps/release/6.3/emboss/apps/transeq.html
    """

    def __init__(self):
        super().__init__('emboss_transeq')

    def __call__(self, seq, frame='1', trim=False, **kwargs):
        """Translate a nucleotide sequence.

        Args:
            seq (str): Sequence to be translated.
                Transeq allows various formats among which FASTA and EMBL
                are most relevant here.
            frame (str): which frames to translate. Defaults to 6, meaning all.
                Other options are: 1, 2, 3, -1, -2, -3, F (all forward) and R.
            trim (bool): Whether transeq should remove trailing * or X.

        Returns:
             WSResponse

        Raises:
            NoDataError: If translation fails or any errors using webservice.
        """
        data = {'sequence': seq,
                'frame': frame,
                'trim': trim,
                }
        for k, v in kwargs.items():
            data[k] = v
        return self._run(data)


class FASTASimilaritySearch(EBIWebService):
    """Use FASTA suite to search a protein sequence against a database.

    The suite contains the different programs - fasta, fastx, fasty, ssearch,
    ggsearch, glsearch, tfastx, tfasty.

    fasta: Scan a protein or DNA sequence library for similar sequences.
        protein:protein or DNA:DNA
    fastx: Compares a DNA sequence to a protein sequence database, translating
        the DNA sequence in three forward (or reverse) frames and allowing
        frameshifts.
        DNA:protein
    fasty: Compares a DNA sequence to a protein sequence database, translating
        the DNA sequence in three forward (or reverse) frames and allowing
        frameshifts.
        DNA:protein
    ssearch: Compare a protein or DNA sequence to a sequence database using
        the Smith-Waterman algorithm.
        protein:protein or DNA:DNA
    ggsearch: Compare a protein or DNA sequence to a sequence database using
        a global alignment (Needleman-Wunsch).
        protein:protein or DNA:DNA
    glsearch: Compare a protein or DNA sequence to a sequence database with
        alignments that are global in the query and local in the database
        sequence (global-local).
        protein:protein or DNA:DNA
    tfastx: Compares a protein sequence to a DNA sequence or DNA sequence
        library. The DNA sequence is translated in three forward and three
        reverse frames, and the protein query sequence is compared to each
        of the six derived protein sequences. The DNA sequence is translated
        from one end to the other; no attempt is made to edit out intervening
        sequences. Termination codons are translated into unknown ('X') amino
        acids.
    """

    def __init__(self):
        super().__init__('fasta')

    def __call__(self, seq, program='fasta', database='uniprotkb_swissprot', stype='protein', **kwargs):
        """Run a sequence similarity search using FASTA suite tools.

        The FASTA tools take many more arguments than listed below which are only the
        mandatory ones. Any further keyword arguments passed to the __call__() method
        will be added to the parameters of the FASTA run.

        Args:
            seq (str): Query sequence(s); preferably in FASTA format.
            database (str): Target database to search against. Default: uniprotkb_swissprot.
                Many choices like uniprotkb, uniprotkb_swissprot, uniprotkb_trembl. For a
                full list run get_target_databases().
            stype (str): Type of query sequence. Default: protein.
                Can be: protein, dna, rna.

        Returns:
            WSResponse
        """
        data = {'sequence': seq,
                'program': program,
                'database': database,
                'stype': stype,
                }
        for k, v in kwargs.items():
            data[k] = v
        return self._run(data)


class SeqChecksum(EBIWebService):
    """Generate checksums for protein/nucleotide sequences."""

    def __init__(self):
        super().__init__('seqcksum')

    def __call__(self, seq, cksmethod='spcrc', stype='protein', **kwargs):
        """Generate checksums for protein/nucleotide sequences.

        The default method generates a CRC64 checksum like it is used in
        UniProtKB.

        Args:
            seq (str): Sequence in FASTA format.
            cksmethod (str): Checksum method. Defaults to: CRC64-ISO (spcrc).
                There are others, including other flavours of CRC, which
                I am not sure about.
            stype (str): Sequence type, protein or dna. Default: protein.

        Returns:
            WSResponse
        """
        data = {'sequence': seq,
                'cksmethod': cksmethod,
                'stype': stype,
                }
        data.update(kwargs)
        return self._run(data)



fasta_search = FASTASimilaritySearch()
translate = EmbossTranseq()
checksum = SeqChecksum()