import hashlib
import re
from collections import defaultdict


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

    def __init__(self, val, typ):
        self.value = val
        self.type = typ
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
        for comment in entry.comments:
            typ, value = comment.split(': ')
            body_and_ev = value.split(' {')
            try:
                body, ev = body_and_ev
            except ValueError:
                print('Weird splitting pattern for comment: {} {}'.format(typ, value))
                continue
            else:
                # handle evidences
                evidences = []
                for token in ev.rstrip('}.').split(', '):
                    try:
                        code, source = token.split('|')
                    except ValueError:
                        code, source = token, None
                    evidences.append(Evidence(code=code, source=source))
                # handle statements
                stmts = re.split('\. ', body)
                for stmt in stmts:
                    text = re.split('\(PubMed:', stmt, 1)[0]
                    for ev in evidences:
                        if ev.source in stmt:
                            anno = Annotation(entry.primary_accession,
                                              Statement(text, typ),
                                              evidence=ev)
                            self.add(anno)

    def size(self):
        """Return length of ACollection list."""
        return self.__len__()

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


if __name__ == '__main__':
    from biocuration import uniprot as up
    with open('C:/Users/kpichler/Documents/Python/evidences/entry.txt', 'r', encoding='ascii') as infile:
        entry = list(up.parse_txt_compatible(infile))[0]
        p = APile()
        p.consume(entry)
        for annotation in p:
            print(annotation)

