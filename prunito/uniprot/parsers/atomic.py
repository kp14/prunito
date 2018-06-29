"""Module for representing UniProtKB data as single annotations.

Many of the fields used in UniProtKB are free text and/or multi-labeled
with many evidence tags attached to one data point. Sometimes it would
be beneficial to be able to treat all data points as separate annotations,
i.e., a data point reported by two different papers would be two separate
annotations. I call these atomic annotations.
"""

import hashlib
import re
from collections import defaultdict
from ...utils import UNIPROT_EVIDENCE_REGEX, PUBMED_REGEX


ECO_PTRN = re.compile('ECO:[0-9]{7}')

class Entity(object):
    """Biological entity which has annotations."""

    def __init__(self):
        self.identifier = None
        self.name = None
        self.taxonomy_id = None
        self.organism = None
        self.sequence = None

    @classmethod
    def from_prunito(cls, record):
        """Create instance from a prunito Record object.
        """
        instance = cls()
        instance.identifier = record.primary_accession
        instance.name = record.recommended_full_name
        instance.organism = record.organism
        instance.taxonomy_id = record.taxonomy_id
        instance.sequence = record.sequence
        return instance

    def __str__(self):
        return '{0}|{1}'.format(self.identifier, self.taxonomy_id)


class Evidence(object):
    """Evidence as defined by Evidence Code Ontology.
    
    Args:
        code (str, optional): Evidence code from ontology; defaults to ECO:0000000.
        source (optional): Source of the evidence, usually a PMID.
        id: md5 signature
    """

    def __init__(self, code='ECO:0000000', source=None):
        if not re.match(ECO_PTRN, code):
            self.code='ECO:0000000'
            # raise ValueError('Invalid ECO code: {}'.format(code))
        else:
            self.code = code
        self.source = source
        self.id = self._calculate_id()

    def is_sim(self):
        return self.code == 'ECO:0000250'

    def is_exp(self):
        return self.code == 'ECO:0000269'

    def _is_auto(self):
        return self.code in ['ECO:0000256', 'ECO:0000259', 'ECO:0000244']

    def __str__(self):
        return '{}-{}'.format(str(self.code), str(self.source))

    def _calculate_id(self):
        """Calculate ID, a digest of object values.

        Returns:
            id, string
        """
        s = self.__str__()
        return hashlib.md5(s.encode()).hexdigest()


class Statement(object):
    """A statement or assertion that can be made about an entity.
    
    Args:
        val (str): The statement.
        typ (str): Type (category) for the statement; taken from UniProt.
    
    Attributes:
        val (str): The statement.
        typ (str): Type (category) for the statement; taken from UniProt.
        id: md5 signature
        """

    def __init__(self, val, typ, **kwargs):
        self.value = val
        self.type = typ
        for k, v in kwargs.items():
            self.__setattr__(k, v)
        self.id = self._calculate_id()

    def __str__(self):
        return '{}-{}'.format(self.type, self.value)

    def _calculate_id(self):
        """Calculate digest of value and type.
        
        Returns:
            id, string
        """
        s = self.__str__()
        return hashlib.md5(s.encode()).hexdigest()


class Annotation(object):
    """An annotation as used in UniProtKB.
    
    Args:
        entity (:obj:`Entity`): The entity an annotation is attached to, a UniProtKB accession.
        stmnt (:obj:`Statement`): The value/type of the annotation statement, e.g. Has xyz activity.
        evidence (:obj:`Evidence`, optional): Evidence for the statement.
    
    Attributes:
        entity (str): The entity an annotation is attached to, a UniProtKB accession.
        evidence (:obj:`Evidence`, optional): Evidence for the statement.
        """

    def __init__(self, entity, stmnt, evidence=None):
        self.entity = entity
        self._statement = stmnt
        self.evidence = evidence
        self.id = self._calculate_id()

    @property
    def value(self):
        """Return value of Annotation statement."""
        return self._statement.value

    @property
    def type(self):
        """Return type of Annotation statement."""
        return self._statement.type

    @property
    def source(self):
        """Return source of Evidence for Annotation statement."""
        try:
            return self.evidence.source
        except AttributeError:
            return None

    @property
    def evidence_code(self):
        """Return the ECO evidence code of an Annotation."""
        try:
            return self.evidence.code
        except AttributeError:
            return None

    def __str__(self):
        return '{en}: {st} - {ev}'.format(en=self.entity.__str__(),
                                          st=self._statement.__str__(),
                                          ev=self.evidence.__str__(),
                                          )

    def __eq__(self, other):
        return self._statement.id == other._statement.id

    def _calculate_id(self):
        """Calculate digest of value and type.

        Returns:
            id, string
        """
        s = self.__str__()
        return hashlib.md5(s.encode()).hexdigest()


class APile(object):
    """A collection (pile) of annotations."""

    def __init__(self):
        self._annotations = []
        self.keywords = defaultdict(list)
        self.rp_tokens = defaultdict(list)

    @classmethod
    def from_iterable(cls, iterable):
        """Alternative constructor to generate an APile.
        
        This assumes an iterable of Annotations instances.
        
        Args:
            iterable: Iterable of Annotation instances.
            
        Returns:
            ACollection instance
        """
        instance = cls()
        for a in iterable:
            instance.add(a)
        return instance

    def add(self, annotation):
        """Add an Annotation to the ACollection`s list."""
        self._annotations.append(annotation)

    def consume(self, entry):
        """Convert a Biopython-type record into Annotations.

        Args:
            entry: Biopython-type Record instance of a UniProt entry
        """
        ap = AtomicParser(entry)
        for annotation in ap.parse():
            self.add(annotation)
        for kw in entry.keywords:
            self.keywords[kw].append(entry.primary_accession)

    def size(self):
        """Return length of ACollection list."""
        return self.__len__()

    def entities(self):
        """Set of entities for which APile contains annotations.

        Returns:
            set
        """
        entities = set()
        for anno in self._annotations:
            entities.add(anno.entity)
        return sorted(entities)

    def sources(self):
        """Set of sources used in evidence tags contained in annotations.

        Returns:
            set
        """
        sources = set()
        for anno in self._annotations:
            if anno.source:
                sources.add(anno.source)
        return sorted(sources)

    def evtags(self):
        """Set of evidence tags contained in annotations.

        Returns:
            set
        """
        evtags = set()
        for anno in self._annotations:
            if anno.evidence_code:
                evtags.add(anno.evidence_code)
        return sorted(evtags)

    def annotation_types(self):
        """Set of evidence tags contained in annotations.

        Returns:
            set
        """
        atypes = set()
        for anno in self._annotations:
            atypes.add(anno.type)
        return sorted(atypes)

    def get_idx(self, idx):
        """Return Annotation at index idx.
        
        Args:
            idx (int): 
        
        Returns:
            Annotation instance
        """
        return self._annotations[idx]

    def __len__(self):
        return self._annotations.__len__()

    def __iter__(self):
        return iter(self._annotations)


class AtomicParser():
    """Parse Biopython-type record entries into atomic annotations.

    Args:
        entry: Biopython-type Record instance of a UniProt entry
    Returns:
        Annotation instances
    """

    def __init__(self, entry):
        self.entry = entry
        self.entity = Entity.from_prunito(entry)
        self._mapper = {'intera': self._parse_interaction,
                       'subcel': self._parse_subcellular_location,
                       'cofact': self._parse_cofactor,
                       'domain': self._parse_freetext,
                       'ptm': self._parse_freetext,
                       }

    def parse(self):
        for comment in self.entry.comments:
            typ, value = comment.split(': ', maxsplit=1)
            parser_func = self._mapper.get(typ[:6].lower(), self._parse_freetext)
            try:
                for annotation in parser_func(typ, value):
                    yield annotation
            except TypeError as e:
                print(e, typ, value)
        for feature in self.entry.features:
            for annotation in self._parse_feature(feature):
                yield annotation
        for kw in self.entity.keywords:
            yield Annotation(self.entity,
                             Statement(kw, 'keyword'),
                             )

    def _parse_feature(self, feature):
        typ, start, stop, description, _ = feature
        evidences = extract_evidences(description)
        text = extract_ft_description(description)
        value = typ + ' ' + text
        if not evidences:
            yield Annotation(self.entity,
                             Statement(value, typ))
        else:
            for e in evidences:
                yield Annotation(self.entity,
                                 Statement(value, typ),
                                 evidence=e)

    def _parse_freetext(self, typ, value):
        """Extract Annotations from freetext comments.

        Args:
            typ (str): type of UniProt comment
            value (str): freetext body of a UniProt comment

        Return:
              Annotation instances
        """
        evidences = extract_evidences(value)
        # There might not be any tags, we provide a token one
        if not evidences:
            evidences.append(Evidence())
        # handling body of annotations
        body = value.split(' {')[0]
        statements = re.split('\. ', body)
        from_pubmed = []
        from_sim = []
        from_unclear = []
        for statement in statements:
                if 'PubMed:' in statement:
                    from_pubmed.append(statement)
                elif '(By similarity)' in statement:
                    from_sim.append(statement)
                elif '(Probable)' in statement:
                    from_pubmed.append(statement)
                else:
                    from_unclear.append(statement)
        for s in from_pubmed:
            text = s.split(' (PubMed:')[0]
            for ev in evidences:
                if ev.source in s:
                    yield Annotation(self.entity,
                                      Statement(text, typ),
                                      evidence=ev)
        for s in from_sim:
            text = s.split(' (By sim')[0]
            sim_tags = [tag for tag in evidences if tag.is_sim()]
            if len(sim_tags) == 1:
                yield Annotation(self.entity,
                                 Statement(text, typ),
                                 evidence=sim_tags[0])
            else:
                yield Annotation(self.entity,
                                 Statement(text, typ),
                                 evidence=Evidence(code='ECO:0000250'))
        for s in from_unclear:
            text = s.rstrip('. ')
            if not from_sim and not from_pubmed and len(evidences) == 1:# 1 tag applicable to all
                yield Annotation(self.entity,
                                 Statement(text, typ),
                                 evidence=evidences[0])
            elif len(evidences) > 1 and len(set([e.code for e in evidences])) == 1:# 1 type of tag applicable to all but unclear source
                this_code = evidences[0].code
                yield Annotation(self.entity,
                                 Statement(text, typ),
                                 evidence=Evidence(code=this_code))
            elif len(evidences) > 1 and len(set([e.code for e in evidences])) > 1 and from_sim:
                code = list(set([e.code for e in evidences if not e.code == 'ECO:0000250']))[0]# a hack
                yield Annotation(self.entity,
                                 Statement(text, typ),
                                 evidence=Evidence(code=code))
            else:
                yield Annotation(self.entity,
                                 Statement(text, typ))

    def _automatic_tags(self, list_of_evidences):
        auto_only = [tag for tag in list_of_evidences if tag.is_auto()]
        return auto_only

    def _sim_tags(self, list_of_evidences):
        sim_only = [tag for tag in list_of_evidences if tag.is_sim()]
        return len(sim_only)

    def _exp_or_opi_tags(self, list_of_evidences):
        exp_only = [tag for tag in list_of_evidences if tag.is_exp()]
        return len(exp_only)

    def _find_tag_for_pubmed(self, pubmed, list_of_evidences ):
        for tag in list_of_evidences:
            if tag.source == pubmed:
                return tag
                break

    def _parse_subcellular_location(self, typ, value):
        """Extract Annotations from subcellular location comments.

        Args:
            typ (str): type of UniProt comment
            value (str): body of a UniProt comment

        Return:
              Annotation instances
        """
        note = '' #TODO: Handle the Note
        try:
            locations_all, note = value.split('. Note=')
        except ValueError: # there is no Note
            locations = value.split('. ')
        else:
            locations = locations_all.split('. ')
        if note:
            yield from self._parse_freetext('SUBCELLULAR LOCATION', note)
        for location in locations:
            try:
                loc, evs = location.split(' {')
            except ValueError:
                loc = location
                evs = None
            if evs:
                evidences = []
                for token in evs.rstrip('}.').split(', '):
                    try:
                        code, source = token.split('|')
                    except ValueError:
                        code, source = token, None
                    evidences.append(Evidence(code=code, source=source))
                for ev in evidences:
                    yield Annotation(self.entity,
                                     Statement(loc, typ),
                                     evidence=ev)
            else:
                yield Annotation(self.entity,
                                 Statement(loc, typ))


    def _parse_cofactor(self, typ, value):
        """Extract Annotations from cofactor comments.

        Args:
            typ (str): type of UniProt comment
            value (str): body of a UniProt comment

        Return:
              Annotation instances
        """
        pass


    def _parse_interaction(self, typ, value):
        """Extract Annotations from interaction comments.

            Args:
                typ (str): type of UniProt comment
                value (str): body of a UniProt comment

            Return:
                  Annotation instances
            """
        pass

    def _parse_keywords(self):
        pass


def extract_evidences(blob):
    evidences = []
    for match in re.finditer(UNIPROT_EVIDENCE_REGEX, blob):
        ev_string = match.group()
        try:
            code, source = ev_string.split('|')
        except ValueError:
            code, source = ev_string, None
        evidences.append(Evidence(code=code, source=source))
    return evidences


def extract_ft_description(blob):
    if blob.startswith('{'):
        return ''
    else:
        try:
            text, _ = blob.split('. {')
        except ValueError:
            return blob
        else:
            return text

