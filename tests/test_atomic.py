import hashlib
import pytest
from prunito.uniprot import atomic


def test_actual_in_expected(actual, expected_list):
    assert actual in expected_list


def test_expected_in_actual(expected, actual_list):
    assert expected in actual_list

# EVS = [('ECO:0000269', '1234567'),
#        ('ECO:0000303', '2345678'),
#        ('ECO:0000305', '33456'),
#        ('ECO:0000269', 'Ref.2'),
#        ]
#
#
# DIGEST = hashlib.md5('{}-{}'.format(EVS[0][0], EVS[0][1]).encode()).hexdigest()
#
#
# VALS = ['Some text here.',
#         'Another sentence.',
#         'A little bit longer this time and no trailing full stop',
#         'A fourth statement.',
#         ]
#
#
# TYPES = ['FUNCTION',
#          'DOUBLE WORD',
#          'MASS SPECTROMETRY',
#          'SUBUNIT',
#          ]
#
#
# ENTITIES = ['A12345',
#             'B12345',
#             'C12345',
#             'D12345',
#             ]
#
#
# def generate_annotations():
#     annotations = []
#     selectors = [(0, 0, 0, 0),
#                  (0, 1, 0, 0),
#                  (0, 2, 1, 2),
#                  (1, 0, 1, 0),
#                  (1, 0, 2, 0),
#                  (2, 1, 0, 1),
#                  (2, 0, 0, 2),
#                  (3, 2, 1, 0),
#                  (3, 0, 2, 1),]
#     for s in selectors:
#         ent = ENTITIES[s[0]]
#         s_val = VALS[s[1]]
#         s_typ = TYPES[s[2]]
#         ev = EVS[s[3]]
#         statement = atomic.Statement(s_val, s_typ)
#         evidence = atomic.Evidence(code=ev[0], source=ev[1])
#         anno = atomic.Annotation(ent, statement, evidence=evidence)
#         annotations.append(anno)
#     return annotations
#
#
# @pytest.fixture
# def ev():
#     return atomic.Evidence(code=EVS[0][0], source=EVS[0][1])
#
#
# @pytest.fixture
# def stmnt():
#     return atomic.Statement(VALS[0], TYPES[0])
#
#
# @pytest.fixture
# def anno():
#     st = atomic.Statement(VALS[1], TYPES[1])
#     ev = atomic.Evidence(code=EVS[0][0], source=EVS[0][1])
#     return atomic.Annotation(ENTITIES[0], st, evidence=ev)
#
#
# @pytest.fixture
# def acoll():
#     data = generate_annotations()
#     coll = atomic.APile.from_iterable(data)
#     return coll
#
#
# def test_ev_code_invalid_eco():
#     with pytest.raises(ValueError):
#         ev = atomic.Evidence(code='ECO:269')
#
#
# def test_ev_id_generation(ev):
#     assert ev.id == DIGEST
#
#
# def test_ev_return_code(ev):
#     assert ev.code == EVS[0][0]
#
#
# def test_ev_return_source(ev):
#     assert ev.source == EVS[0][1]
#
#
# def test_st_id_swapped_input(stmnt):
#     s = '{}-{}'.format(VALS[0], TYPES[0])
#     hash = hashlib.md5(s.encode()).hexdigest()
#     assert stmnt.id != hash
#
#
# def test_st_id(stmnt):
#     s = '{}-{}'.format(TYPES[0], VALS[0])
#     hash = hashlib.md5(s.encode()).hexdigest()
#     assert stmnt.id == hash
#
#
# def test_anno_value(anno):
#     assert anno.value == VALS[1]
#
#
# def test_anno_source(anno):
#     assert anno.source == EVS[0][1]
#
#
# def test_anno_type(anno):
#     assert anno.type == TYPES[1]
#
#
# def test_anno_entity(anno):
#     assert anno.entity == ENTITIES[0]
#
#
# def test_anno_ev_code(anno):
#     assert anno.evidence_code == 'ECO:0000269'
#
#
# def test_anno_str_representation(anno):
#     expected = '{en}: {stt}-{stv} - {evc}-{evs}'.format(en=ENTITIES[0],
#                                                         stt=TYPES[1],
#                                                         stv=VALS[1],
#                                                         evc=EVS[0][0],
#                                                         evs=EVS[0][1])
#     assert anno.__str__() == expected
#
#
# def test_collection_constructor_alternative(acoll):
#     assert acoll.size() == 9
#
#
# def test_anno_compare(acoll):
#     a1 = acoll.get_idx(0)
#     a2 =acoll.get_idx(6)
#     assert a1 == a2
