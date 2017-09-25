import hashlib
import re
from biocuration.utils import UNIPROT_EVIDENCE_REGEX, PUBMED_REGEX


ECO_PTRN = re.compile('ECO:[0-9]{7}')

class Evidence(object):
    """Evidence as defined by Evidence Code Ontology.
    
    Args:
        code (str, optional): Evidence code from ontology; defaults to ECO:0000000.
        source (optional): Source of the evidence, usually a PMID.
        id: md5 signature
    """

    def __init__(self, code='ECO:0000000', source=None):
        if not re.match(ECO_PTRN, code):
            raise ValueError('Invalid ECO code: {}'.format(code))
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
        entity (str): The entity an annotation is attached to, a UniProtKB accession.
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
        return '{en}: {st} - {ev}'.format(en=self.entity,
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
        return entities

    def sources(self):
        """Set of sources used in evidence tags contained in annotations.

        Returns:
            set
        """
        sources = set()
        for anno in self._annotations:
            sources.add(anno.source)
        return sources

    def evtags(self):
        """Set of evidence tags contained in annotations.

        Returns:
            set
        """
        evtags = set()
        for anno in self._annotations:
            evtags.add(anno.evidence_code)
        return evtags

    def annotation_types(self):
        """Set of evidence tags contained in annotations.

        Returns:
            set
        """
        atypes = set()
        for anno in self._annotations:
            atypes.add(anno.type)
        return atypes

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

    def _parse_feature(self, feature):
        typ, start, stop, description, _ = feature
        try:
            text, evs = description.split('. {')
        except ValueError:
            text = description
            evs = None
        value = typ + ' ' + text
        if evs:
            evidences = []
            for token in evs.rstrip('}.').split(', '):
                try:
                    code, source = token.split('|')
                except ValueError:
                    code, source = token, None
                evidences.append(Evidence(code=code, source=source))
            for ev in evidences:
                yield Annotation(self.entry.primary_accession,
                                 Statement(value, typ),
                                 evidence=ev)
        else:
            yield Annotation(self.entry.primary_accession,
                             Statement(value, typ))

    def _parse_freetext(self, typ, value):
        """Extract Annotations from freetext comments.

        Args:
            typ (str): type of UniProt comment
            value (str): freetext body of a UniProt comment

        Return:
              Annotation instances
        """
        evidences = self._extract_evidences(value)
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
                    yield Annotation(self.entry.primary_accession,
                                      Statement(text, typ),
                                      evidence=ev)
        for s in from_sim:
            text = s.split(' (By sim')[0]
            sim_tags = [tag for tag in evidences if tag.is_sim()]
            if len(sim_tags) == 1:
                yield Annotation(self.entry.primary_accession,
                                 Statement(text, typ),
                                 evidence=sim_tags[0])
            else:
                yield Annotation(self.entry.primary_accession,
                                 Statement(text, typ),
                                 evidence=Evidence(code='ECO:0000250'))
        for s in from_unclear:
            pass



    def _extract_evidences(self, blob):
        evidences = []
        for match in re.finditer(UNIPROT_EVIDENCE_REGEX, blob):
            ev_string = match.group()
            try:
                code, source = ev_string.split('|')
            except ValueError:
                code, source = ev_string, None
            evidences.append(Evidence(code=code, source=source))
        return evidences

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

        # body_and_ev = value.split(' {')
        # try:
        #     body, ev = body_and_ev
        # except ValueError:
        #     print('Weird splitting pattern for comment: {} {}'.format(typ, value))
        # else:
        #     # handle evidences
        #     evidences = []
        #     for token in ev.rstrip('}.').split(', '):
        #         try:
        #             code, source = token.split('|')
        #         except ValueError:
        #             code, source = token, None
        #         evidences.append(Evidence(code=code, source=source))
        #     # handle statements
        #     stmts = re.split('\. ', body)
        #     for stmt in stmts:
        #         text = re.split('\(PubMed:', stmt, 1)[0]
        #         if len(evidences) == 1:
        #             anno = Annotation(self.entry.primary_accession,
        #                               Statement(text, typ),
        #                               evidence=ev)
        #             yield anno
        #         else:
        #             for ev in evidences:
        #                 if ev.source in stmt:
        #                     anno = Annotation(self.entry.primary_accession,
        #                                       Statement(text, typ),
        #                                       evidence=ev)
        #                     yield anno


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
                    yield Annotation(self.entry.primary_accession,
                                     Statement(loc, typ),
                                     evidence=ev)
            else:
                yield Annotation(self.entry.primary_accession,
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

if __name__ == '__main__':
    from biocuration import uniprot as up
    with open('C:/users/kpichler/Documents/Python/biocuration_pc/tests/SwissProt/atomic.txl', 'r', encoding='ascii') as infile:
    # with open('/home/klemens/Documents/data.txt', 'r', encoding='ascii') as infile:
    # with open('/home/klemens/Downloads/entry.txt', 'r', encoding='ascii') as infile:
        # entry = list(up.parse_txt_compatible(infile))[0]
        p = APile()
        for entry in up.parse_txt_compatible(infile):
            p.consume(entry)
        for annotation in p:
            print(annotation)


