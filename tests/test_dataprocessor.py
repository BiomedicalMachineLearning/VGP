#!/usr/bin/env python

"""Tests for `stlearn` package."""


import unittest

from transprs.datasets.classes import DataProcessor
import pandas as pd
import numpy as np


class TestDataProcessor(unittest.TestCase):
    """Tests for Preprocessing step."""

    def test_DataProcessor(self):

        sumstats = pd.read_csv("tests/test_data/test_sumstats.csv")
        imputed_snps = np.array(pd.read_csv("tests/test_data/test_bim.csv")["1"])

        processor = DataProcessor(sumstats=sumstats, imputed_snps=imputed_snps)

        processor.clean_snps()
        processor.filter_imputed(info=0.9)
        processor.extract_intersection()
        processor.split_chromosomes()
        processor.sort_snps_chr()
