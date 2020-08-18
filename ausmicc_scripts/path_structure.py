#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and functions related to path structure

Created on Aug 18 2020

@author: V.R.Marcelino
"""
import os

# path to ref databases:
class ref_db():
    def __init__(self,blastdb_16Smicro="/home/ref_databases/BlastDBs/16SMicrobial/16SMicrobial"):

        self.blastdb_16Smicro = blastdb_16Smicro


# folder structure:
class amplicon_paths():
    def __init__(self, abi_collec="/mnt/datastore/AMPLICON/ABI",
                 fasta_collec="/mnt/datastore/AMPLICON/FASTA"):
        
        self.abi_collec = abi_collec
        self.fasta_collec = fasta_collec


## make folders to store AUSMICC genomes and amplicons
def make_db_folders():
        os.mkdir("/mnt/datastore")
        os.mkdir("/mnt/datastore/AMPLICON")
        os.mkdir("/mnt/datastore/GENOMES")
        os.mkdir("/mnt/datastore/AMPLICON/ABI")
        os.mkdir("/mnt/datastore/AMPLICON/FASTA")
        os.mkdir("/mnt/datastore/GENOMES/RAW")
        os.mkdir("/mnt/datastore/GENOMES/PAIRED")
        os.mkdir("/mnt/datastore/GENOMES/CONTIGS")
        os.mkdir("/mnt/datastore/GENOMES/ANNOTATIONS")


