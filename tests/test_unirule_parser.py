# -*- coding: utf-8 -*-
"""
Tests for the UniRule parser.

Created on Thu Feb  4 15:06:48 2016

@author: kp14
"""

import pytest
import os

from prunito.uniprot import parse_rules

# Some diffs are very long here if tests fail so we set this
# to None to display the whole diff
# TestCase.maxDiff = None

# Import paths for data files
base_dir = os.path.dirname(__file__)
filename = "selected-rules.xml"
data_file = os.path.join(base_dir, "UniRule", filename)

rule = list(parse_rules(data_file))[0]


def test_meta_id():
    assert rule.meta["id"] == "UR000172772"


def test_meta_creator():
    assert rule.meta["creator"] == "amcdouga"


def test_meta_mod_by():
    assert rule.meta["modifiedBy"] == "h.zellner"


def test_meta_status():
    assert rule.meta["status"] == "APPLY"


def test_meta_created():
    assert rule.meta["created"] == "2015-08-03+01:00"


def test_meta_old_id():
    assert rule.meta["oldRuleNum"] == "RU362117"


def test_meta_uniprot_id():
    assert rule.meta["uniprotId"] == "Cyt_b"


def test_meta_dataclass():
    assert rule.meta["dataClass"] == "Protein"


def test_length_conditions_sets_main():
    assert len(rule.main.conditions) == 1


def test_main_cond_pfam():
    assert rule.main.conditions[0][0].type == "Pfam id"
    assert rule.main.conditions[0][1].type == "Pfam id"


def test_main_cond_iter_negative():
    for c in rule.iter_conditions():
        assert c.negative == False


def test_main_cond_taxon():
    assert rule.main.conditions[0][2].type == "taxon"
    assert rule.main.conditions[0][2].value == "Eukaryota"
    assert rule.main.conditions[0][2].cvId == "2759"


def test_main_cond_mito():
    assert rule.main.conditions[0][3].type == "gene location"
    assert rule.main.conditions[0][3].value == "Mitochondrion"
    assert rule.main.conditions[0][3].negative == False


def test_number_annotations_main():
    assert len(rule.main.annotations) == 13


def test_type_annotations_main():
    assert rule.main.annotations[0].type == "recommendedName"


def test_class_annotations_main():
    assert rule.main.annotations[0].class_ == "protein"


def test_number_of_cases():
    assert len(rule.cases) == 6


def test_all_case_conditions_values():
    expected = [
        "Fungi",
        "Vertebrata",
        "Fungi",
        "Saccharomyces",
        "Saccharomyces",
        "Fungi",
        "Saccharomyces",
    ]
    parsed = [c.value for c in rule.iter_case_conditions()]
    assert parsed == expected


def test_all_case_conditions_negativity():
    expected = [False, False, True, False, True, False, True]
    parsed = [c.negative for c in rule.iter_case_conditions()]
    assert parsed == expected


def test_case_1_annotation():
    parsed = rule.cases[0].annotations[0]
    assert parsed.value == "TRM5"
    assert parsed.type == "primary"
    assert parsed.class_ == "gene"


def test_case_2_annotation():
    parsed = rule.cases[1].annotations[0]
    assert parsed.value == "TRMT5"
    assert parsed.type == "primary"
    assert parsed.class_ == "gene"
    parsed = rule.cases[1].annotations[1]
    assert parsed.value == "TRM5"
    assert parsed.type == "synonym"
    assert parsed.class_ == "gene"


def test_case_3_annotation():
    for a in rule.cases[2].annotations:
        assert a.class_ == "protein"
        assert a.type in ["recommendedName", "alternativeName"]
        assert a.subtype in ["fullName", "ecNumber"]


def test_sam_numbers():
    assert len(rule.sam_features) == 3


def test_sam_triggers():
    map_ = {0: "transmembrane", 1: "signal", 2: "coiledCoil"}
    for idx, sam in enumerate(rule.sam_features):
        assert sam.trigger == map_[idx]


def test_sam_0_str():
    expected = (
        "Number of conditions: 1\n"
        "Number of annotations: 0\n"
        "Trigger: transmembrane\n"
        "Range: 1-12\n"
    )
    assert rule.sam_features[0].__str__() == expected
