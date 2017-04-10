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

    def _calculate_id(self):
        """Calculate ID, a digest of object values.

        Returns:
            id, string
        """
        s = '{}-{}'.format(str(self.code), str(self.source))
        return hashlib.md5(s.encode()).hexdigest()


class Statement(object):
    """A statement or assertion that can be made about an entity."""

    def __init__(self, val, typ):
        self.value = val
        self.type = typ
        self.id = self._calculate_id()

    def _calculate_id(self):
        """Calculate digest of value and type.
        
        Returns:
            id, string
        """
        s = '{}-{}'.format(self.type, self.value)
        return hashlib.md5(s.encode()).hexdigest()


class Annotation(object):
    """An annotation is used in UniProtKB."""

    def __init__(self, stmnt, entity, evidence=None):
        self._statement = stmnt
        self.entity = entity
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
        self._statements = {}
        self._state_map = defaultdict(set)

    def add(self, annotation):
        self._statements[annotation.digest] = annotation._statement
        annotation._statement = annotation.digest
        self._state_map[annotation].add(annotation)


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
