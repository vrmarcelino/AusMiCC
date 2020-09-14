#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and functions related to path structure

Created on Aug 18 2020

@author: V.R.Marcelino
"""
from pathlib import Path

# path to reference databases:
class ref_db():
    def __init__(self,blastdb_16Smicro="/home/ref_databases/BlastDBs/16SMicrobial/16SMicrobial",
                 blastdb_silva="/home/ref_databases/BlastDBs/16S_SILVA_132_NCBI.fna"):

        self.blastdb_16Smicro = blastdb_16Smicro
        self.blastdb_silva = blastdb_silva


# folder structure:
class amplicon_paths():
    def __init__(self, ab1_collec="/mnt/datastore/AMPLICON/AB1",
                 fasta_collec="/mnt/datastore/AMPLICON/FASTA"):
        
        self.ab1_collec = ab1_collec
        self.fasta_collec = fasta_collec


class genome_paths():
    def __init__(self, raw_reads="/mnt/datastore/GENOMES/RAW",
                 paired="/mnt/datastore/GENOMES/PAIRED",
                 unpaired="/mnt/datastore/GENOMES/UNPAIRED",
                 contigs="/mnt/datastore/GENOMES/CONTIGS",
                 annotations="/mnt/datastore/GENOMES/ANNOTATIONS"):
        self.raw_reads = raw_reads
        self.paired = paired
        self.unpaired = unpaired
        self.contigs = contigs
        self.annotations = annotations


# make folders to store AUSMICC genomes and amplicons
# similar names to what is found in mba
def make_db_folders():
        Path("/mnt/datastore/AMPLICON").mkdir(parents=True, exist_ok=True)
        Path("/mnt/datastore/AMPLICON/AB1").mkdir(parents=True, exist_ok=True)
        Path("/mnt/datastore/AMPLICON/FASTA").mkdir(parents=True, exist_ok=True)
        Path("/mnt/datastore/GENOMES/RAW").mkdir(parents=True, exist_ok=True)
        Path("/mnt/datastore/GENOMES/PAIRED").mkdir(parents=True, exist_ok=True)
        Path("/mnt/datastore/GENOMES/UNPAIRED").mkdir(parents=True, exist_ok=True)
        Path("/mnt/datastore/GENOMES/CONTIGS").mkdir(parents=True, exist_ok=True)
        Path("/mnt/datastore/GENOMES/ANNOTATIONS").mkdir(parents=True, exist_ok=True)
        print ("\n Inexistent folders created. \n")


