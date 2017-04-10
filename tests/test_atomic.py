import hashlib
import pytest
from biocuration.uniprot import atomic

CODE = 'ECO:0000269'
SOURCE = '1234567'
EVS = [('ECO:0000269', '1234567'),
       ('ECO:0000269', '2345678'),
       ('ECO:0000269', '33456'),
       ('ECO:0000269', 'Ref.2'),
       ]

DIGEST = hashlib.md5('{}-{}'.format(CODE, SOURCE).encode()).hexdigest()

VALS = ['Some text here.',
        'Another sentence.',
        'A little bit longer this time and no trailing full stop',
        ]

TYPES = ['FUNCTION',
         'DOUBLE WORD',
         'MASS SPECTROMETRY',
         ]

ENTITIES = ['A12345',
            'B12345',
            'C12345',
            ]

@pytest.fixture
def ev():
    return atomic.Evidence(code=EVS[0][0], source=EVS[0][1])

@pytest.fixture
def stmnt():
    return atomic.Statement(VALS[0], TYPES[0])

@pytest.fixture
def anno():
    st = atomic.Statement(VALS[1], TYPES[1])
    ev = atomic.Evidence(code=EVS[0][0], source=EVS[0][1])
    return atomic.Annotation(ENTITIES[0], st, evidence=ev)

def test_ev_code_invalid_eco():
    with pytest.raises(ValueError):
        ev = atomic.Evidence(code='ECO:269')

def test_ev_id_generation(ev):
    assert ev.id == DIGEST

def test_ev_return_code(ev):
    assert ev.code == EVS[0][0]

def test_ev_return_source(ev):
    assert ev.source == EVS[0][1]

def test_st_id_swapped_input(stmnt):
    s = '{}-{}'.format(VALS[0], TYPES[0])
    hash = hashlib.md5(s.encode()).hexdigest()
    assert stmnt.id != hash

def test_st_id(stmnt):
    s = '{}-{}'.format(TYPES[0], VALS[0])
    hash = hashlib.md5(s.encode()).hexdigest()
    assert stmnt.id == hash

def test_anno_value(anno):
    assert anno.value == VALS[1]

def test_anno_source(anno):
    assert anno.source == SOURCE

def test_anno_type(anno):
    assert anno.type == TYPES[1]

def test_anno_entity(anno):
    assert anno.entity == ENTITIES[0]