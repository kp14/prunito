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
filename = 'made_up.xml'
data_file = os.path.join(base_dir, 'UniRule', filename)

rule = list(parse_rules(data_file))[0]

def test_meta_info():
    assert(rule.meta['id'], 'UR000150961')
    assert(rule.meta['creator'], 'kpichler')
    assert(rule.meta['status'], 'TEST')


def test_length_conditions_sets_main():
    assert(len(rule.main.conditions), 3)


def test_number_annotations_main():
    assert(len(rule.main.annotations), 2)


def test_type_annotations_main():
    assert(rule.main.annotations[0]['type_'], 'subcellularLocation')


def test_class_annotations_main():
    assert(rule.main.annotations[0]['class_'], 'comment')
