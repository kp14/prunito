# coding: utf-8
#import functools
import itertools
import logging
import re
from collections import defaultdict
from .atomic import APile
from ..utils import EC_REGEX


logging.basicConfig(level=logging.WARN)

LINE_ENDINGS = ' \n\r'


class Record:
    '''Class emulating biopython.SwissProt.Record.

    All fields biopython.SwissProt.Record provides are provided, too.
    Addinitional ones also, for convenience.
    '''
    def __init__(self, bag):
        self._bag = bag
        self._features = []

    @property
    def entry_name(self):
        return self._bag['ID'][0][0]

    @property
    def data_class(self):
        return self._bag['ID'][0][1].strip(';')

    @property
    def molecule_type(self):
        return None

    @property
    def sequence_length(self):
        return int(self._bag['ID'][0][2])

    @property
    def accessions(self):
        return list(itertools.chain(*self._bag['AC']))

    @property
    def created(self):
        date, _ = self._bag['DT'][0].strip('.').split(', ')
        return (date, 0) # Version is always zero here

    @property
    def sequence_update(self):
        date, version_string = self._bag['DT'][1].strip('.').split(', ')
        version_split = version_string.split()
        return (date, int(version_split[2]))

    @property
    def annotation_update(self):
        date, version_string = self._bag['DT'][2].strip('.').split(', ')
        version_split = version_string.split()
        return (date, int(version_split[2]))

    @property
    def primary_accession(self):
        return self._bag['AC'][0][0]

    @property
    def description(self):
        return ' '.join(self._bag['DE'])

    @property
    def recommended_full_name(self):
        """Full name of protein as annotated in UniProtKB.

        In the TXT format, this corresponds to the value of line:
        DE   RecName: Full= ...;
        Evidence tags are stripped out.
        """
        rec_line =  self._bag['DE'][0].strip(' ;')
        assert rec_line.startswith('RecName')
        name_part = rec_line.split(' {')[0]
        rec_name = name_part[14:]
        return rec_name

    @property
    def ec_numbers(self):
        """Get all EC numbers present in the DE block.

        Returns:
            list of str
        """
        return re.findall(EC_REGEX, self.description)

    @property
    def gene_name(self):

        return ' '.join(self._bag['GN'])

    @property
    def primary_gene_name(self):
        """Gene name as captured in UniProtKB.

        Specifically, the Name= field in GN lines. If an entry does not have
        a Name field even if it has other GN fields this returns an empty string.

        Returns:
            str
        """
        primary_gn = ''
        try:
            first_gn_token = self.gene_name.split(';')[0]
        except KeyError:
            pass
        else:
            if first_gn_token.startswith('Name'):
                primary_gn = first_gn_token[5:].split(' {')[0]
        return primary_gn

    @property
    def organism(self):
        return ' '.join(self._bag['OS'])

    @property
    def organelle(self):
        try:
            return self._bag['OG'][0]
        except IndexError:
            return ''

    @property
    def organism_classification(self):
        return _flatten_lists(self._bag['OC'])

    @property
    def taxonomy_id(self):
        return self._bag['OX']

    @property
    def host_organism(self):
        orgs = []
        for oh_line in self._bag['OH']:
            _, org = oh_line.strip('.').split('; ')
            orgs.append(org)
        return orgs

    @property
    def host_taxonomy_id(self):
        ids = []
        for oh_line in self._bag['OH']:
            id_string, _ = oh_line.strip('.').split('; ')
            _, taxid = id_string.split('=')
            ids.append(taxid)
        return ids

    @property
    def references(self):
        ref_list = []
        for num, data in self._bag['RF'].items():
            ref_list.append(Reference(num, data))
        return ref_list

    @property
    def all_pubmed_ids(self):
        """Return all PMIDs found in the references of an entry."""
        id_list = []
        for ref in self.references:
          for source in ref.references:
              name, identifier = source
              if name == 'PubMed':
                  id_list.append(identifier)
        return id_list

    @property
    #@functools.lru_cache(maxsize=1)
    def comments(self):
        concat = ' '.join(self._bag['CC'])
        topics = concat.split('-!- ')
        topics = [x.strip() for x in topics]
        return topics[1:]

    @property
    def cross_references(self):
        return self._bag['DR']

    @property
    def keywords(self):
        return _flatten_lists(self._bag['KW'])

    @property
    #@functools.lru_cache(maxsize=1)
    def seqinfo(self):
        tokens = self._bag['SQ'][0].split()
        length = int(tokens[1])
        mol_weight = int(tokens[3])
        crc64 = tokens[5]
        return (length, mol_weight, crc64)

    @property
    #@functools.lru_cache(maxsize=1)
    def features(self):
        if self._features:
            return self._features
        else:
            current = None
            ft_bag_len = len(self._bag['FT'])
            for idx, featureline in enumerate(self._bag['FT']):
                logging.debug('Looking at featureline: {}'.format(featureline))
                if featureline[0]:
                    if current:
                        self._features.append(tuple(current))
                        current = featureline.copy()
                    else:
                        current  = featureline.copy()
                else:
                    if featureline[3].startswith('/FTId'):
                        current[4] = featureline[3][6:-1]
                    else:
                        if current[3].endswith('-'):
                            add_space = ''
                        else:
                            add_space = ' '
                        current[3] += '{0}{1}'.format(add_space, featureline[3])
                if idx == ft_bag_len - 1:
                    self._features.append(tuple(current))
            return self._features

    @property
    def sequence(self):
        gapped_seq = ''.join(self._bag['  '])
        seq = gapped_seq.replace(' ', '')
        return seq

    @property
    def isoforms(self):
        """
        Generates isoforms.
        :return: tuple (name, sequence)
        """
        slen = len(self.sequence)
        iso = defaultdict(dict)
        for ft in self.features:
            if ft[0] == 'VAR_SEQ':
                for match in re.findall('isoform [A-Za-z0-9-]+', ft[3]):
                    name = match[8:]
                    try:
                        _ = iso[name]['blocks']
                    except KeyError:
                        iso[name]['blocks'] = []
                    if ft[3].startswith('Missing'):
                        iso[name]['blocks'].append((ft[1] - 1, ft[2], '-'))
                    else:
                        _, change = ft[3].split(' -> ')
                        change_no_iso = re.sub(' [(]in isoform .+[])]\.', '', change)
                        change_no_tags = re.sub(' [{].+[}]\.', '', change_no_iso)
                        aa_changes = change_no_tags.replace(' ', '')
                        iso[name]['blocks'].append((ft[1] - 1, ft[2], aa_changes))
        for isoform, data in iso.items():
            seq = [0]
            for block in data['blocks']:
                seq.append(block[0])
                seq.append(block[1])
            seq.append(slen)
            for idx in range(0, len(seq), 2):
                start = seq[idx]
                stop = seq[idx + 1]
                if start >= stop:
                    pass
                else:
                    tpl = (start, stop, '+')
                    data['blocks'].append(tpl)
            data['blocks'].sort()
            temp = []
            for block in data['blocks']:
                if block[2] == '+':
                    for idx in range(block[0], block[1]):
                        temp.append(self.sequence[idx])
                elif block[2] == '-':
                    pass
                else:
                    temp.extend(list(block[2]))
            data['seq'] = ''.join(temp)
        for iso_name, data in iso.items():
            qualified_iso_name = '{0}-{1}'.format(self.primary_accession, iso_name)
            yield (qualified_iso_name, data['seq'])


class Reference():

    def __init__(self, num, data):
        self._num = num
        self._data = data

    @property
    def number(self):
        return self._num

    @property
    def positions(self):
        return self._data['RP']
        # pos = ' '.join(self._data['RP'])
        # return pos.strip('.')

    @property
    def comments(self):
        return self._parse_rc_rx('RC')

    @property
    def references(self):
        return self._parse_rc_rx('RX')

    @property
    def pmid(self):
        """
        Return PubMed ID of reference, None if there is none.
        :return: PubMed ID, string
        """
        for source, identifier in self.references:
            if source == 'PubMed':
                return identifier

    @property
    def authors(self):
        if len(self._data['RA']) == 1:
            authors = self._data['RA'][0]
        else:
            authors = ' '.join(self._data['RA'])
        return authors.rstrip(';')

    @property
    def first_author(self):
        """
        Return first author name of reference.
        :return: author name, string
        """
        authors = self.authors.split(', ')
        return authors[0]

    @property
    def title(self):
        title = ' '.join(self._data['RT'])
        return title.strip('";')

    @property
    def location(self):
        return self._data['RL'][0]

    def _parse_rc_rx(self, context):
        context = context
        if len(self._data[context]) == 1:
            comments = self._data[context][0]
        else:
            comments = ' '.join(self._data[context])
        tokens = comments.split('; ')
        c_tuples = []
        for tok in tokens:
            if tok:
                key, val = tok.rstrip(';').split('=')
                c_tuples.append((key, val))
        return c_tuples

def _flatten_lists(lol):
    '''Flatten list of lists.'''
    flat = []
    for sublist in lol:
        flat.extend(sublist)
    return flat


def parse_txt(handle):
    bag, context = _set_up()
    if handle:
        ignore = False
        for line in handle:
            line_type, line = line[:2], line[5:]
            if line_type in ['id', 'dt', 'ac']:
                ignore = True
            stripped_line = line.rstrip()
    #        stripped_line = line
            try:
                if line_type.startswith('R'):
                    if line_type == 'RN':
                        number = _extract_ref_number(line)
                        context = number
                        bag['RF'][context] = defaultdict(list)
                    else:
                        bag['RF'][context][line_type].append(stripped_line)
                value = PARSER_MAP[line_type](stripped_line)
                if value:
                    bag[line_type].append(value)
            except KeyError:
                if line_type == '//':
                    if not ignore:
                        yield bag
                        bag, context = _set_up()
                    else:
                        logging.warning('Ignored lowercase entry.')
                        del bag
                        bag, context = _set_up()
                        ignore = False
                else:
                    logging.debug("Unknown line type: {}".format(line_type))
    else:
        print('Nothing to parse. Are you sure the query returns anything?')


def parse_txt_compatible(handle):
    for bag in parse_txt(handle):
        yield Record(bag)


def parse_txt_atomic(handle):
    pile = APile()
    for entry in parse_txt_compatible(handle):
        pile.consume(entry)
    return pile

def _set_up():
    bag = defaultdict(list)
    bag['RF'] = {}
    context = None
    return bag, context

def _extract_ref_number(line):
    pattern = re.compile('(\[[0-9]+\])')
    tokens = re.split(pattern, line)
    number = int(tokens[1].strip('[]'))
    return number

def _parse_generic(line):
    return line


def _parse_id(line):
    id_line_tokens = line.split()
    return id_line_tokens


def _parse_accessions(line):
    stripped_line = line.strip(';.')
    tokens = stripped_line.split(';')
    accs = [x.strip() for x in tokens]
    return accs


def _parse_cc(line):
    return line.strip()

def _parse_species(line):
    """Extract species information from OS line

    :line: TODO
    :returns: TODO

    """
    species = line.split(' (')[0]
    return species


def _parse_taxonomy(line):
    stripped_line = line.strip(';.')
    tokens = stripped_line.split(';')
    nodes = [x.strip() for x in tokens]
    return nodes


def _parse_taxid(line):
    """Extract NCBI taxonomy ID.

    :line: TODO
    :returns: TODO

    """
    _, id_string = line.strip(';').split('=')
    taxid = id_string.split()[0]
    return taxid


def _parse_cc(line):
    stripped_line = line.strip()
    if not stripped_line[0:3] in ['---', 'Cop', 'Dis']:
        return stripped_line


def _parse_dr(line):
    stripped_line = line.strip(';.')
    tokens = stripped_line.split('; ')
    return tuple(tokens)


def _parse_kw(line):
    stripped_line = line.strip(';.')
    tokens = stripped_line.split(';')
    kwds = [x.strip() for x in tokens]
    return kwds


def _parse_ft(line):
    '''Parses components of a UniProt flatfile FT line.

    FT lines follow the following format:
    [FT   ]KEY      x      y       Description.
    Assuming a maximum value for x or y of 99999, which seems OK as
    the longest sequences are about 35000 amino acids. X and y can
    be numbers or of format >123 etc.

    Returns:
    list
    '''
    key = line[0:8].strip()
    start_ft = _try_parsing_as_int(line[10:15].strip())
    end_ft = _try_parsing_as_int(line[17:22].strip())
    description = line[29:].strip()
    # if not key:
    #     if description.startswith('/FT'):
    #         _parse_ft.previous[-1] = description[6:-1]
    #     else:
    #         _parse_ft.previous[3] += new_desc
    featureline = [key, start_ft, end_ft, description, '']
    logging.debug('Parsed FT: {}'.format(featureline))
    return featureline


def _try_parsing_as_int(token):
    try:
        return int(token.strip())
    except ValueError:
        return token.strip()


PARSER_MAP = {"ID": _parse_id,
              "AC": _parse_accessions,
              "DT": _parse_generic,
              "DE": _parse_cc,
              "GN": _parse_generic,
              "OS": _parse_generic,
              "OG": _parse_generic,
              "OC": _parse_taxonomy,
              "OH": _parse_generic,
              "OX": _parse_taxid,
              "RF": _parse_generic,
              "CC": _parse_cc,
              "DR": _parse_dr,
              "KW": _parse_kw,
              "PE": _parse_generic,
              "FT": _parse_ft,
              "**": _parse_generic,
              "SQ": _parse_generic,
              "  ": _parse_generic}

if __name__ == '__main__':
    datafile = 'C:/Users/kpichler/Documents/Python/biocuration_pc/tests/SwissProt/one_sp_entry.txl'
    with open(datafile, 'r', encoding='ascii') as data:
        e = list(parse_txt_compatible(data))
        e[0].features