import hashlib
import pytest
from biocuration.uniprot import atomic

CODE = 'ECO:0000269'
SOURCE = '1234567'
DIGEST = hashlib.md5('{}-{}'.format(CODE, SOURCE).encode()).hexdigest()
VALS = ['Some text here.',
        'Another sentence.',
        'A little bit longer this time and no trailing full stop',
        ]
TYPES = ['FUNCTION',
         'DOUBLE WORD',
         'MASS SPECTROMETRY',
         ]

@pytest.fixture
def ev():
    return atomic.Evidence(code=CODE, source=SOURCE)

@pytest.fixture
def stmnt():
    return atomic.Statement(VALS[0], TYPES[0])

def test_ev_code_invalid_eco():
    with pytest.raises(ValueError):
        ev = atomic.Evidence(code='ECO:269')

def test_ev_id_generation(ev):
    assert ev.id == DIGEST

def test_ev_return_code(ev):
    assert ev.code == CODE

def test_ev_return_source(ev):
    assert ev.source == SOURCE

def test_st_id_swapped_input(stmnt):
    s = '{}-{}'.format(VALS[0], TYPES[0])
    hash = hashlib.md5(s.encode()).hexdigest()
    assert stmnt.id != hash

def test_st_id(stmnt):
    s = '{}-{}'.format(TYPES[0], VALS[0])
    hash = hashlib.md5(s.encode()).hexdigest()
    assert stmnt.id == hash