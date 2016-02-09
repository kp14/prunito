# -*- coding: utf-8 -*-
"""
Tests for the UniRule parser.

Created on Thu Feb  4 15:06:48 2016

@author: kp14
"""

import pytest
import os

from biocuration.uniprot import parse_rules

# Some diffs are very long here if tests fail so we set this
# to None to display the whole diff
#TestCase.maxDiff = None

# Import paths for data files
base_dir = os.path.dirname(__file__)
filename = 'selected-rules.xml'
data_file = os.path.join(base_dir, 'UniRule', filename)

rule = list(parse_rules(data_file))[0]

def test_meta_id():
    assert(rule.meta['id'], 'UR000172772')


def test_meta_creator():
    assert(rule.meta['creator'], 'amcdouga')

def test_meta_mod_by():
    assert(rule.meta['modifiedBy'], 'h.zellner')

def test_meta_status():
    assert(rule.meta['status'], 'APPLY')

def test_meta_created():
    assert(rule.meta['created'], '2015-08-13+01:00')

def test_meta_old_id():
    assert(rule.meta['oldRuleNum'], 'RU362117')

def test_meta_uniprot_id():
    assert(rule.meta['uniprotId'], 'Cyt_b')

def test_meta_dataclass():
    assert(rule.meta['dataClass'], 'Protein')


def test_length_conditions_sets_main():
    assert(len(rule.main.conditions), 1)

def test_main_cond_pfam():
    assert(rule.main.conditions[0][0].type, 'Pfam id')
    assert(rule.main.conditions[0][1].type, 'Pfam id')

def test_main_cond_iter_negative():
    for c in rule.iter_conditions():
        assert(c.negative, False)

def test_main_cond_taxon():
    assert(rule.main.conditions[0][2].type, 'taxon')
    assert(rule.main.conditions[0][2].value, 'Eukaryota')

def test_main_cond_mito():
    assert(rule.main.conditions[0][3].type, 'gene_location')
    assert(rule.main.conditions[0][3].value, 'Mitochondrion')
    assert(rule.main.conditions[0][3].negative, False)


def test_number_annotations_main():
    assert(len(rule.main.annotations), 13)


def test_type_annotations_main():
    assert(rule.main.annotations[0].type, 'recommendedName')


def test_class_annotations_main():
    assert(rule.main.annotations[0].class_, 'protein')
