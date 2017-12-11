from prunito import uniprot as up


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


# def test_tax_ids_info_several_ids():
#     expected = {
#               "taxonomies": [
#                 {
#                   "taxonomyId": 9606,
#                   "mnemonic": "HUMAN",
#                   "scientificName": "Homo sapiens",
#                   "commonName": "Human",
#                   "rank": "species",
#                   "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/9605",
#                   "childrenLinks": [
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/741158",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/63221"
#                   ],
#                   "siblingsLinks": [
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1425170"
#                   ]
#                 },
#                 {
#                   "taxonomyId": 3202,
#                   "mnemonic": "9MARC",
#                   "scientificName": "Jungermannia",
#                   "rank": "genus",
#                   "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/3201",
#                   "childrenLinks": [
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/402631",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/37392",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/402630",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1670782",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1670783",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1670784",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/350782",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/588650",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/588651",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1112838",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/362817",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/362816",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/463576",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/362815",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/362814",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/248332",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/746483",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/280831",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/746484",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/746485",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/280830",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/746486",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/3203"
#                   ],
#                   "siblingsLinks": [
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/402090",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/53014",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/209812",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1331044",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1867305",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/984532",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/306424",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/209808",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1527784",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/463570",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/248327"
#                   ]
#                 },
#                 {
#                   "taxonomyId": 176652,
#                   "mnemonic": "IIV6",
#                   "scientificName": "Invertebrate iridescent virus 6",
#                   "commonName": "IIV-6",
#                   "synonym": "Chilo iridescent virus",
#                   "rank": "species",
#                   "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10487",
#                   "childrenLinks": [
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/132417"
#                   ],
#                   "siblingsLinks": [
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/327984",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/176651"
#                   ]
#                 },
#                 {
#                   "taxonomyId": 10090,
#                   "mnemonic": "MOUSE",
#                   "scientificName": "Mus musculus",
#                   "commonName": "Mouse",
#                   "rank": "species",
#                   "parentLink": "https://www.ebi.ac.uk/proteins/api/taxonomy/id/862507",
#                   "childrenLinks": [
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/477815",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/477816",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/179238",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/35531",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1879032",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10091",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10092",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/57486",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/39442",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/947985",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1266728",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1643390",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/80274",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/46456",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/116058",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/1385377"
#                   ],
#                   "siblingsLinks": [
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10096",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10097",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10098",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10100",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10089",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/473865",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/254704",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/186842",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/10103",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/186193",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/27681",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/481680",
#                     "https://www.ebi.ac.uk/proteins/api/taxonomy/id/83773"
#                   ]
#                 }
#               ]
#             }
#
#     actual = up.tax_node_info('9606,3202,176652,10090')
#     assert len(expected['taxonomies']) == len(actual['taxonomies'])
#     for tax in zip(expected['taxonomies'], actual['taxonomies']):
#         for k, v in tax[0].items():
#             tax[1][k] == v