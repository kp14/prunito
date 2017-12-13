import pytest
from prunito import uniprot as up
from prunito.uniprot.rest.proteinsapi import NoDataError


def test_get_info_on_taxID_string_id():
    expected = {
              "taxonomyId": 9606,
              "mnemonic": "HUMAN",
              "scientificName": "Homo sapiens",
              "commonName": "Human",
              "rank": "species",
              "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/9605",
              "childrenLinks": [
                "https://www.ebi.ac.uk/proteins/api/taxonomy/id/741158",
                "https://www.ebi.ac.uk/proteins/api/taxonomy/id/63221"
              ],
              "siblingsLinks": [
                "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1425170"
              ]
            }
    actual = up.get_info_on_taxID('9606') # taxID passed as str
    for k, v in expected.items():
        assert actual[k] == v


def test_get_info_on_taxID_invalid_taxid():
    with pytest.raises(NoDataError):
        r = up.get_info_on_taxID('*89HJ')


def test_get_info_on_taxID_nonexisting_taxid():
    with pytest.raises(NoDataError):
        r = up.get_info_on_taxID(10000001000000)


def test_get_info_on_taxID_int_id():
    expected = {
              "taxonomyId": 9606,
              "mnemonic": "HUMAN",
              "scientificName": "Homo sapiens",
              "commonName": "Human",
              "rank": "species",
              "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/9605",
              "childrenLinks": [
                "https://www.ebi.ac.uk/proteins/api/taxonomy/id/741158",
                "https://www.ebi.ac.uk/proteins/api/taxonomy/id/63221"
              ],
              "siblingsLinks": [
                "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1425170"
              ]
            }
    actual = up.get_info_on_taxID(9606) # taxID passed as int
    for k, v in expected.items():
        assert actual[k] == v


def test_tax_ids_info_several_ids_list():
    expected = [
                {
                  "taxonomyId": 9606,
                  "mnemonic": "HUMAN",
                  "scientificName": "Homo sapiens",
                  "commonName": "Human",
                  "rank": "species",
                  "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/9605",
                  "childrenLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/741158",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/63221"
                  ],
                  "siblingsLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1425170"
                  ]
                },
                {
                  "taxonomyId": 3202,
                  "mnemonic": "9MARC",
                  "scientificName": "Jungermannia",
                  "rank": "genus",
                  "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/3201",
                  "childrenLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/402631",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/37392",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/402630",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1670782",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1670783",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1670784",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/350782",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/588650",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/588651",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1112838",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/362817",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/362816",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/463576",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/362815",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/362814",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/248332",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/746483",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/280831",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/746484",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/746485",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/280830",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/746486",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/3203"
                  ],
                  "siblingsLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/402090",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/53014",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/209812",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1331044",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1867305",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/984532",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/306424",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/209808",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1527784",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/463570",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/248327"
                  ]
                },
                {
                  "taxonomyId": 176652,
                  "mnemonic": "IIV6",
                  "scientificName": "Invertebrate iridescent virus 6",
                  "commonName": "IIV-6",
                  "synonym": "Chilo iridescent virus",
                  "rank": "species",
                  "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10487",
                  "childrenLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/132417"
                  ],
                  "siblingsLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/327984",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/176651"
                  ]
                },
                {
                  "taxonomyId": 10090,
                  "mnemonic": "MOUSE",
                  "scientificName": "Mus musculus",
                  "commonName": "Mouse",
                  "rank": "species",
                  "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/862507",
                  "childrenLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/477815",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/477816",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/179238",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/35531",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1879032",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10091",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10092",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/57486",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/39442",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/947985",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1266728",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1643390",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/80274",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/46456",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/116058",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1385377"
                  ],
                  "siblingsLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10096",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10097",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10098",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10100",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10089",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/473865",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/254704",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/186842",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10103",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/186193",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/27681",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/481680",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/83773"
                  ]
                }
              ]

    actual = up.get_info_on_taxIDs([9606,3202,176652,10090])
    for tax in zip(expected, list(actual)):
        for k, v in tax[0].items():
            assert tax[1][k] == v


def test_tax_ids_info_several_ids_tuple():
    expected = [
                {
                  "taxonomyId": 9606,
                  "mnemonic": "HUMAN",
                  "scientificName": "Homo sapiens",
                  "commonName": "Human",
                  "rank": "species",
                  "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/9605",
                  "childrenLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/741158",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/63221"
                  ],
                  "siblingsLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1425170"
                  ]
                },
                {
                  "taxonomyId": 3202,
                  "mnemonic": "9MARC",
                  "scientificName": "Jungermannia",
                  "rank": "genus",
                  "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/3201",
                  "childrenLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/402631",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/37392",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/402630",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1670782",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1670783",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1670784",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/350782",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/588650",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/588651",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1112838",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/362817",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/362816",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/463576",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/362815",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/362814",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/248332",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/746483",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/280831",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/746484",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/746485",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/280830",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/746486",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/3203"
                  ],
                  "siblingsLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/402090",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/53014",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/209812",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1331044",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1867305",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/984532",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/306424",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/209808",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1527784",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/463570",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/248327"
                  ]
                },
                {
                  "taxonomyId": 176652,
                  "mnemonic": "IIV6",
                  "scientificName": "Invertebrate iridescent virus 6",
                  "commonName": "IIV-6",
                  "synonym": "Chilo iridescent virus",
                  "rank": "species",
                  "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10487",
                  "childrenLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/132417"
                  ],
                  "siblingsLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/327984",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/176651"
                  ]
                },
                {
                  "taxonomyId": 10090,
                  "mnemonic": "MOUSE",
                  "scientificName": "Mus musculus",
                  "commonName": "Mouse",
                  "rank": "species",
                  "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/862507",
                  "childrenLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/477815",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/477816",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/179238",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/35531",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1879032",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10091",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10092",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/57486",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/39442",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/947985",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1266728",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1643390",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/80274",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/46456",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/116058",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1385377"
                  ],
                  "siblingsLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10096",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10097",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10098",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10100",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10089",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/473865",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/254704",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/186842",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10103",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/186193",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/27681",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/481680",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/83773"
                  ]
                }
              ]
    actual = up.get_info_on_taxIDs((9606,3202,176652,10090))
    for tax in zip(expected, list(actual)):
        for k, v in tax[0].items():
            assert tax[1][k] == v


def test_tax_ids_info_several_ids_set():
    expected = [
                {
                  "taxonomyId": 9606,
                  "mnemonic": "HUMAN",
                  "scientificName": "Homo sapiens",
                  "commonName": "Human",
                  "rank": "species",
                  "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/9605",
                  "childrenLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/741158",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/63221"
                  ],
                  "siblingsLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1425170"
                  ]
                },
                {
                  "taxonomyId": 3202,
                  "mnemonic": "9MARC",
                  "scientificName": "Jungermannia",
                  "rank": "genus",
                  "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/3201",
                  "childrenLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/402631",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/37392",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/402630",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1670782",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1670783",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1670784",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/350782",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/588650",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/588651",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1112838",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/362817",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/362816",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/463576",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/362815",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/362814",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/248332",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/746483",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/280831",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/746484",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/746485",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/280830",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/746486",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/3203"
                  ],
                  "siblingsLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/402090",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/53014",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/209812",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1331044",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1867305",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/984532",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/306424",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/209808",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1527784",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/463570",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/248327"
                  ]
                },
                {
                  "taxonomyId": 176652,
                  "mnemonic": "IIV6",
                  "scientificName": "Invertebrate iridescent virus 6",
                  "commonName": "IIV-6",
                  "synonym": "Chilo iridescent virus",
                  "rank": "species",
                  "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10487",
                  "childrenLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/132417"
                  ],
                  "siblingsLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/327984",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/176651"
                  ]
                },
                {
                  "taxonomyId": 10090,
                  "mnemonic": "MOUSE",
                  "scientificName": "Mus musculus",
                  "commonName": "Mouse",
                  "rank": "species",
                  "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/862507",
                  "childrenLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/477815",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/477816",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/179238",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/35531",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1879032",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10091",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10092",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/57486",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/39442",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/947985",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1266728",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1643390",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/80274",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/46456",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/116058",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1385377"
                  ],
                  "siblingsLinks": [
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10096",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10097",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10098",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10100",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10089",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/473865",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/254704",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/186842",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10103",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/186193",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/27681",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/481680",
                    "https://www.ebi.ac.uk/proteins/api/taxonomy/id/83773"
                  ]
                }
              ]
    actual = up.get_info_on_taxIDs({9606,3202,176652,10090})
    for tax in list(actual):
        assert tax['taxonomyId'] in {9606,3202,176652,10090}


def test_get_info_on_taxIDs_string_input():
    r = up.get_info_on_taxIDs('9606,10090')
    with pytest.raises(ValueError):
        next(r)


def test_get_info_on_taxIDs_invalid_taxid():
    r = up.get_info_on_taxIDs(['8UGH5'])
    with pytest.raises(NoDataError):
        next(r)


def test_get_info_on_taxIDs_nonexisting_taxid():
    r = up.get_info_on_taxIDs([960600000])
    with pytest.raises(NoDataError):
        next(r)


def test_get_lineage_for_taxID():
    expected = [
                {
                  "taxonomyId": 9606,
                  "scientificName": "Homo sapiens"
                },
                {
                  "taxonomyId": 9605,
                  "scientificName": "Homo"
                },
                {
                  "taxonomyId": 207598,
                  "scientificName": "Homininae"
                },
                {
                  "taxonomyId": 9604,
                  "scientificName": "Hominidae"
                },
                {
                  "taxonomyId": 314295,
                  "scientificName": "Hominoidea"
                },
                {
                  "taxonomyId": 9526,
                  "scientificName": "Catarrhini"
                },
                {
                  "taxonomyId": 314293,
                  "scientificName": "Simiiformes"
                },
                {
                  "taxonomyId": 376913,
                  "scientificName": "Haplorrhini"
                },
                {
                  "taxonomyId": 9443,
                  "scientificName": "Primates"
                },
                {
                  "taxonomyId": 314146,
                  "scientificName": "Euarchontoglires"
                },
                {
                  "taxonomyId": 1437010,
                  "scientificName": "Boreoeutheria"
                },
                {
                  "taxonomyId": 9347,
                  "scientificName": "Eutheria"
                },
                {
                  "taxonomyId": 32525,
                  "scientificName": "Theria"
                },
                {
                  "taxonomyId": 40674,
                  "scientificName": "Mammalia"
                },
                {
                  "taxonomyId": 32524,
                  "scientificName": "Amniota"
                },
                {
                  "taxonomyId": 32523,
                  "scientificName": "Tetrapoda"
                },
                {
                  "taxonomyId": 1338369,
                  "scientificName": "Dipnotetrapodomorpha"
                },
                {
                  "taxonomyId": 8287,
                  "scientificName": "Sarcopterygii"
                },
                {
                  "taxonomyId": 117571,
                  "scientificName": "Euteleostomi"
                },
                {
                  "taxonomyId": 117570,
                  "scientificName": "Teleostomi"
                },
                {
                  "taxonomyId": 7776,
                  "scientificName": "Gnathostomata"
                },
                {
                  "taxonomyId": 7742,
                  "scientificName": "Vertebrata"
                },
                {
                  "taxonomyId": 89593,
                  "scientificName": "Craniata"
                },
                {
                  "taxonomyId": 7711,
                  "scientificName": "Chordata"
                },
                {
                  "taxonomyId": 33511,
                  "scientificName": "Deuterostomia"
                },
                {
                  "taxonomyId": 33213,
                  "scientificName": "Bilateria"
                },
                {
                  "taxonomyId": 6072,
                  "scientificName": "Eumetazoa"
                },
                {
                  "taxonomyId": 33208,
                  "scientificName": "Metazoa"
                },
                {
                  "taxonomyId": 33154,
                  "scientificName": "Opisthokonta"
                },
                {
                  "taxonomyId": 2759,
                  "scientificName": "Eukaryota"
                },
                {
                  "taxonomyId": 131567,
                  "scientificName": "cellular organisms"
                },
                {
                  "taxonomyId": 1,
                  "scientificName": "root"
                }
              ]
    actual = up.get_lineage_for_taxID(9606)
    actual_list = list(actual)
    assert len(expected) == len(actual_list)
    for item in zip(expected, actual_list):
        for k, v in item[0].items():
            assert item[1][k] == v


def test_get_lineage_for_taxID_invalid_taxid():
    r = up.get_lineage_for_taxID('*89HJ')
    with pytest.raises(NoDataError):
        next(r)


def test_get_lineage_for_taxID_nonexisting_taxid():
    r = up.get_lineage_for_taxID(1000010000000)
    with pytest.raises(NoDataError):
        next(r)