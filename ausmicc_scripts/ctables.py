#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class to store info for the different MySQL tables.

Created on Sep 15 2020

@author: V.R.Marcelino
"""

# store info obtained from 16S amplicon sequencing:

class info_16s():
    def __init__(self, isolate_name=None, ab1_loc=None, sequence=None,primer=None,
                 full_len_sequence=None,lca_taxid=None,sp_taxid=None, species=None):

        # goes to 16S table:
        self.isolate_name = isolate_name
        self.ab1_loc = ab1_loc
        self.sequence = sequence
        self.primer = primer

        # goes to isolate table:
        self.full_len_sequence = full_len_sequence
        self.lca_taxid = lca_taxid
        self.sp_taxid = sp_taxid
        self.species = species

