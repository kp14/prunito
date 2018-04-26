#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test for the SwissProt parser on SwissProt files.

Created on Mon May 18 11:59:36 2015

@author: kpichler
"""
import pytest
import os

from Bio import SwissProt
import prunito.uniprot as up


# Import paths for data files
base_dir = os.path.dirname(__file__)
filename = 'one_sp_entry.txl'
#filename = 'many_sp_entries.txl'
#datafile = os.path.join('SwissProt', filename)
datafile = os.path.join(base_dir, 'SwissProt', filename)

with open(datafile, "r", encoding="ascii") as data:
    my_record = list(up.parse_txt(data))


with open(datafile, "r", encoding="ascii") as data:
    biopython_record = SwissProt.read(data)


def test_wrong_format_paramater_for_search():
    with pytest.raises(ValueError):
        up.search('test', frmt='xsd')


def test_uppercase_format_parameter_works_for_search():
    pass


def test_object_is_record_instance():
    assert isinstance(my_record[0], up.parsers.parser_knowledgebase_txt.Record) == True


def test_data_class():
    assert my_record[0].data_class == 'Reviewed'


def test_sequence_length():
    assert my_record[0].sequence_length == 1013


def test_entry_name():
    assert my_record[0].entry_name == "NOD2_BOVIN"


def test_primary_accession():
    assert my_record[0].primary_accession == "Q6E804"


def test_description():
    assert my_record[0].description == "RecName: Full=Nucleotide-binding oligomerization domain-containing protein 2; EC=2.3.1.234 {ECO:0000255|HAMAP-Rule:MF_01445, ECO:0000269|PubMed:22378793}; AltName: Full=Caspase recruitment domain-containing protein 15; EC=6.3.1.999;"


def test_recommended_full_name():
    assert my_record[0].recommended_full_name == 'Nucleotide-binding oligomerization domain-containing protein 2'


def test_ec_numbers():
    assert my_record[0].ec_numbers == ['2.3.1.234', '6.3.1.999']

def test_gene_name():
    assert my_record[0].gene_name == 'Name=NOD2 {ECO:0000255|HAMAP-Rule:MF_01445}; Synonyms=CARD15;'


def test_primary_gene_name():
    assert my_record[0].primary_gene_name == 'NOD2'


def test_organism_classification():
    assert my_record[0].organism_classification == ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Euteleostomi', 'Mammalia', 'Eutheria', 'Laurasiatheria', 'Cetartiodactyla', 'Ruminantia', 'Pecora', 'Bovidae', 'Bovinae', 'Bos']


def test_organism():
    assert my_record[0].organism == 'Bos taurus (Bovine).'


def test_reference_number():
    assert my_record[0].references[0].number == 1


def test_reference_positions():
    assert my_record[0].references[0].positions == ['NUCLEOTIDE SEQUENCE [GENOMIC DNA / MRNA], AND VARIANTS ALA-70;',
                                                         'MET-196; ASN-505; GLN-681; HIS-689; ARG-733 AND LEU-1007.']


def test_reference_comments():
    assert my_record[0].references[0].comments == []


def test_reference_references():
    assert my_record[0].references[0].references == [('PubMed', '16897345'),
                                                     ('DOI', '10.1007/s00335-005-0148-2')]


def test_reference_pmid_retrieval():
    assert my_record[0].references[0].pmid == '16897345'


def test_reference_authors():
    assert my_record[0].references[0].authors == 'Taylor K.H., Taylor J.F., White S.N., Womack J.E.'


def test_reference_first_author():
    assert my_record[0].references[0].first_author == 'Taylor K.H.'


def test_reference_title():
    assert my_record[0].references[0].title == 'Identification of genetic variation and putative regulatory regions in bovine CARD15.'


def test_reference_location():
    assert my_record[0].references[0].location == 'Mamm. Genome 17:892-901(2006).'


def test_list_of_pmids():
    assert my_record[0].all_pubmed_ids == ['16897345', '20698950', '16203728', '18511561', '19592251', '21887730']


def test_comments_number_of():
    assert len(my_record[0].comments) == 9


def test_comments_content():
    assert my_record[0].comments[6] == 'SIMILARITY: Contains 2 CARD domains. {ECO:0000255|PROSITE-ProRule:PRU00046}.'


def test_cross_references():
    assert my_record[0].cross_references[0] == ('EMBL', 'AY518737', 'AAS09824.1', '-', 'mRNA')
    assert my_record[0].cross_references[18] == ('PaxDb', 'Q6E804', '-')
    assert my_record[0].cross_references[26] == ('HOGENOM', 'HOG000113814', '-')


def test_keywords():
    assert my_record[0].keywords[0] == 'ATP-binding'
    assert my_record[0].keywords == ['ATP-binding', 'Cell membrane', 'Complete proteome', 'Cytoplasm', 'Immunity', 'Innate immunity', 'Leucine-rich repeat', 'Membrane', 'Nucleotide-binding', 'Polymorphism', 'Reference proteome', 'Repeat', 'Ubl conjugation']
    assert len(my_record[0].keywords) == 13


def test_seqinfo():
    assert my_record[0].seqinfo == (1013, 112800, 'ACAE5F8C4BEFE11B')


def test_features():
    assert len(my_record[0].features[0]) == 5
    assert my_record[0].features[0] == ('CHAIN', 1, 1013, 'Nucleotide-binding oligomerization domain-containing protein 2.', 'PRO_0000375804')


def test_sequence():
    assert len(my_record[0].sequence) == 1013


def test_entry_name_compatibility():
    assert biopython_record.entry_name == my_record[0].entry_name


def test_accessions_compatibility():
    assert biopython_record.accessions == my_record[0].accessions


def test_organism_classification_compatibility():
    assert biopython_record.organism_classification == my_record[0].organism_classification


def test_organism_compatibility():
    assert biopython_record.organism == my_record[0].organism


def test_sequence_compatibility():
    assert my_record[0].sequence == biopython_record.sequence


def test_description_compatibility():
    assert my_record[0].description == biopython_record.description


def test_gene_name_compatbility():
    assert my_record[0].gene_name == biopython_record.gene_name


def test_primary_accession_compatibility():
    assert my_record[0].primary_accession == biopython_record.accessions[0]


def test_ref_author_compatibility():
    assert my_record[0].references[0].authors == biopython_record.references[0].authors


def test_comment_compatibility():
    assert len(biopython_record.comments) == len(my_record[0].comments)
    for i in range(0, len(biopython_record.comments)):
        assert biopython_record.comments[i] == my_record[0].comments[i]


def test_feature_compatibility():
    assert biopython_record.features[0] == my_record[0].features[0]
    assert len(biopython_record.features) == len(my_record[0].features)


def test_ignore_lowercase_entries():
    datafile = os.path.join(base_dir, 'SwissProt', 'contains_lowercase.txl')
    with open(datafile, 'r') as infile:
        entries = list(up.parse_txt(infile))
        assert len(entries) == 23