import hashlib
import pytest
from biocuration.uniprot import atomic


EXPECTED_ANNOTATIONS = [
    'P62140: FUNCTION-Dephosphorylates the 'Ser-418' residue of FOXP3 in regulatory T- cells (Treg) from patients with rheumatoid arthritis, thereby inactivating FOXP3 and rendering Treg cells functionally defective - ECO:0000269-PubMed:23396208',
    'P62140: FUNCTION-Protein phosphatase that associates with over 200 regulatory proteins to form highly specific holoenzymes which dephosphorylate hundreds of biological targets - ECO:0000269-None',
    'P62140: FUNCTION-Protein phosphatase (PP1) is essential for cell division, it participates in the regulation of glycogen metabolism, muscle contractility and protein synthesis - ECO:0000269-None',
    'P62140: FUNCTION-Involved in regulation of ionic conductances and long-term synaptic plasticity - ECO:0000269-None',
    'P62140: FUNCTION-Component of the PTW/PP1 phosphatase complex, which plays a role in the control of chromatin structure and cell cycle progression during the transition from mitosis into interphase - ECO:0000269-None',
    'P62140: FUNCTION-In balance with CSNK1D and CSNK1E, determines the circadian period length, through the regulation of the speed and rhythmicity of PER1 and PER2 phosphorylation - ECO:0000269-None',
    'P62140: FUNCTION-May dephosphorylate CSNK1D and CSNK1E  - ECO:0000269-None',
    'P62140: CATALYTIC ACTIVITY-[a protein]-serine/threonine phosphate + H(2)O = [a protein]-serine/threonine + phosphate - None',
    'P62140: CATALYTIC ACTIVITY-[Myosin light-chain] phosphate + H(2)O = [myosin light-chain] + phosphate - ECO:0000269-PubMed:23396208',
    'P62140: SUBUNIT-Interacts with RRP1B  - ECO:0000269-PubMed:20926688',
    'P62140: SUBCELLULAR LOCATION-Cytoplasm - ECO:0000269-PubMed:11739654',
    'P62140: SUBCELLULAR LOCATION-Nucleus - ECO:0000269-PubMed:11739654',
    'P62140: SUBCELLULAR LOCATION-Nucleus, nucleoplasm - ECO:0000269-PubMed:11739654',
    'P62140: SUBCELLULAR LOCATION-Nucleus, nucleolus - ECO:0000269-PubMed:11739654',
    'P62140: SUBCELLULAR LOCATION-Nucleus, nucleolus - ECO:0000269-PubMed:20926688',
    'P62140: SUBCELLULAR LOCATION-Note=Highly mobile in cells and can be relocalized through interaction with targeting subunits. In the presence of PPP1R8 relocalizes from the nucleus to nuclear speckles - None'
    'P62140: INDUCTION-Up-regulated in synovial fluid mononuclear cells and peripheral blood mononuclear cells from patients with rheumatoid arthritis - ECO:0000269-PubMed:23396208',
    'P62140: SIMILARITY-Belongs to the PPP phosphatase family - ECO:0000305-None',
    'P62140: SIMILARITY-PP-1 subfamily - ECO:0000305-None',
    'P62140: INIT_MET-INIT_MET Removed - ECO:0000244-PubMed:22223895',
    'P62140: INIT_MET-INIT_MET Removed - ECO:0000244-PubMed:22814378',
    'P62140: INIT_MET-INIT_MET Removed - ECO:0000269-PubMed:12665801',
    'P62140: CHAIN-CHAIN Serine/threonine-protein phosphatase PP1-beta catalytic subunit. - None',
    'P62140: ACT_SITE-ACT_SITE Proton donor - ECO:0000250-None',
    'P62140: METAL-METAL Manganese 1 - ECO:0000250-None',
    'P62140: METAL-METAL Manganese 1 - ECO:0000250-None',
    'P62140: METAL-METAL Manganese 1 - ECO:0000250-None',
    'P62140: METAL-METAL Manganese 2 - ECO:0000250-None',
    'P62140: METAL-METAL Manganese 2 - ECO:0000250-None',
    'P62140: METAL-METAL Manganese 2 - ECO:0000250-None',
    'P62140: METAL-METAL Manganese 2 - ECO:0000250-None',
    'P62140: MOD_RES-MOD_RES N-acetylalanine - ECO:0000244-PubMed:22223895',
    'P62140: MOD_RES-MOD_RES N-acetylalanine - ECO:0000244-PubMed:22814378',
    'P62140: MOD_RES-MOD_RES N-acetylalanine - ECO:0000269-PubMed:12665801',
    'P62140: MOD_RES-MOD_RES Phosphothreonine - ECO:0000244-PubMed:18669648',
    'P62140: MOD_RES-MOD_RES Phosphothreonine - ECO:0000244-PubMed:23186163',
    'P62140: VARIANT-VARIANT P -> R (probable disease-associated mutation found in patients with rasopathy reminiscent of NSLH) - ECO:0000269-PubMed:27264673',
    'P62140: VARIANT-VARIANT A -> P (probable disease-associated mutation found in patients with rasopathy reminiscent of NSLH) - ECO:0000269-PubMed:27264673',
    'P62140: CONFLICT-CONFLICT L -> P (in Ref. 5; AAV38549) - ECO:0000305-None',
]

EVS = [('ECO:0000269', '1234567'),
       ('ECO:0000303', '2345678'),
       ('ECO:0000305', '33456'),
       ('ECO:0000269', 'Ref.2'),
       ]


DIGEST = hashlib.md5('{}-{}'.format(EVS[0][0], EVS[0][1]).encode()).hexdigest()


VALS = ['Some text here.',
        'Another sentence.',
        'A little bit longer this time and no trailing full stop',
        'A fourth statement.',
        ]


TYPES = ['FUNCTION',
         'DOUBLE WORD',
         'MASS SPECTROMETRY',
         'SUBUNIT',
         ]


ENTITIES = ['A12345',
            'B12345',
            'C12345',
            'D12345',
            ]


def generate_annotations():
    annotations = []
    selectors = [(0, 0, 0, 0),
                 (0, 1, 0, 0),
                 (0, 2, 1, 2),
                 (1, 0, 1, 0),
                 (1, 0, 2, 0),
                 (2, 1, 0, 1),
                 (2, 0, 0, 2),
                 (3, 2, 1, 0),
                 (3, 0, 2, 1),]
    for s in selectors:
        ent = ENTITIES[s[0]]
        s_val = VALS[s[1]]
        s_typ = TYPES[s[2]]
        ev = EVS[s[3]]
        statement = atomic.Statement(s_val, s_typ)
        evidence = atomic.Evidence(code=ev[0], source=ev[1])
        anno = atomic.Annotation(ent, statement, evidence=evidence)
        annotations.append(anno)
    return annotations


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


@pytest.fixture
def acoll():
    data = generate_annotations()
    coll = atomic.APile.from_iterable(data)
    return coll


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
    assert anno.source == EVS[0][1]


def test_anno_type(anno):
    assert anno.type == TYPES[1]


def test_anno_entity(anno):
    assert anno.entity == ENTITIES[0]


def test_anno_ev_code(anno):
    assert anno.evidence_code == 'ECO:0000269'


def test_anno_str_representation(anno):
    expected = '{en}: {stt}-{stv} - {evc}-{evs}'.format(en=ENTITIES[0],
                                                        stt=TYPES[1],
                                                        stv=VALS[1],
                                                        evc=EVS[0][0],
                                                        evs=EVS[0][1])
    assert anno.__str__() == expected


def test_collection_constructor_alternative(acoll):
    assert acoll.size() == 9


def test_anno_compare(acoll):
    a1 = acoll.get_idx(0)
    a2 =acoll.get_idx(6)
    assert a1 == a2
