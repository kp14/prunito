#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test for the SwissProt parser on SwissProt files.

See conftest.py for more details on automated test generation.

Created on Mon May 18 11:59:36 2015

@author: kpichler
"""
from unittest import TestCase

#from Bio import SwissProt
#import biocuration.uniprotkb.parsing as up

# Some diffs are very long here if tests fail so we set this
# to None to display the whole diff
TestCase.maxDiff = None


def test_compare_annotations(my_annotation, biopy_annotation):
    assert my_annotation == biopy_annotation