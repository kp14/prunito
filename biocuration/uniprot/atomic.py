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
    def evidence(self):
        try:
            return self.evidence.code
        except AttributeError:
            return None

    def __str__(self):
        return '{en}: {st} - {ev}'.format(en=self.entity,
                                          st=self._statement.__str__(),
                                          ev=self.evidence.__str__(),
                                          )

    def __cmp__(self, other):
        return self.id == other.id

    def _calculate_id(self):
        """Calculate digest of value and type.

        Returns:
            id, string
        """
        s = '{}-{}-{}'.format(self.entity, self.evidence.id, self._statement.id)
        return hashlib.md5(s.encode()).hexdigest()


class ACollection(object):
    """A collection of annotations."""

    def __init__(self):
        self._annotations = set()

    def add(self, annotation):
        self._annotations.add(annotation)

    def __iter__(self):
        return self

    def __next__(self):
        for anno in self._annotations:
            yield anno



if __name__ == '__main__':
    from biocuration import uniprot as up
    with open('entry.txt', 'r', encoding='ascii') as infile:
        entry = list(up.parse_txt_compatible(infile))[0]
        for c in entry.comments:
            typ, val = c.split(': ')
            body_and_ev = val.split(' {')
            if len(body_and_ev) == 2:
                body, ev = body_and_ev
            else:
                print('Weird splitting pattern for comment: {} {}'.format(typ, val))
            stmts = re.split('\. ', body)
            print(stmts)
            print(ev.rstrip('}.'))
