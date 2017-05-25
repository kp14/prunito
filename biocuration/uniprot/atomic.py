import hashlib
import re
from collections import defaultdict


ECO_PTRN = re.compile('ECO:[0-9]{7}')

class Evidence(object):
    """Evidence as defined by ECO."""

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
    """A statement or assertion that can be made about an entity."""

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
    """An annotation is used in UniProtKB."""

    def __init__(self, entity, stmnt, evidence=None):
        self.entity = entity
        self._statement = stmnt
        self.evidence = evidence
        self.id = self._calculate_id()

    @property
    def value(self):
        return self._statement.value

    @property
    def type(self):
        return self._statement.type

    @property
    def source(self):
        try:
            return self.evidence.source
        except AttributeError:
            return None

    @property
    def evidence_code(self):
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


class ACollection(object):
    """A collection of annotations."""

    def __init__(self):
        self._annotations = []

    @classmethod
    def from_iterable(cls, iterable):
        """Alternative constructor."""
        instance = cls()
        for a in iterable:
            instance.add(a)
        return instance

    def add(self, annotation):
        self._annotations.append(annotation)

    def size(self):
        return self.__len__()

    def get_idx(self,idx):
        """Return Annotation at index idx.
        Returns:
            Annotation instance
        """
        return self._annotations[idx]

    def __len__(self):
        return self._annotations.__len__()

    def __iter__(self):
        return iter(self._annotations)


def consume(entry):
    coll = ACollection()
    for c in entry.comments:
        typ, val = c.split(': ')
        body_and_ev = val.split(' {')
        try:
            body, ev = body_and_ev
        except ValueError:
            print('Weird splitting pattern for comment: {} {}'.format(typ, val))
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
                        coll.add(anno)
    for anno in coll:
        print(anno)


if __name__ == '__main__':
    from biocuration import uniprot as up
    with open('C:/Users/kpichler/Documents/Python/evidences/entry.txt', 'r', encoding='ascii') as infile:
        entry = list(up.parse_txt_compatible(infile))[0]
        consume(entry)

