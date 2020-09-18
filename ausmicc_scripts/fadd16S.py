#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Helper functions to add 16S to the database

Created on 18 Sep  2020

@author: V.R.Marcelino
"""

from ausmicc_scripts import ctables
from ausmicc_scripts import fconnector

def add_16S_record(in_obj, in_cursor):

    # get a dict of isolate_names:ids
    isolateid_map = fconnector.get_ids([in_obj.isolate_name], 'isolate', in_cursor)
    idisolate = isolateid_map[in_obj.isolate_name]

    # convert sequence from string to blob
    partial_sequence = ''.join(format(i, 'b') for i in bytearray(in_obj.sequence, encoding ='utf-8'))

    query = "INSERT INTO info16s (idisolate,ab1_file_loc,sequence,primer) VALUES (%s, %s,%s,%s)"
    val = (idisolate,in_obj.ab1_loc,partial_sequence,in_obj.primer)

    in_cursor.execute(query, val)
    print ("Added %s, %s record to the database\n" %(in_obj.isolate_name, in_obj.primer))

def update_isolate_info(in_obj,in_cursor):
    isolateid_map = fconnector.get_ids([in_obj.isolate_name], 'isolate', in_cursor)
    idisolate = isolateid_map[in_obj.isolate_name]

    #convert full len sequence to bytes (blob), if it exists:
    if in_obj.full_len_sequence != None:
        full_len_sequence = ''.join(format(i, 'b') for i in bytearray(in_obj.full_len_sequence, encoding ='utf-8'))
    else:
        full_len_sequence = None

    query = "UPDATE isolate set full_length_seq=%s,species=%s,species_taxid=%s WHERE idisolate=%s"
    val = (full_len_sequence,in_obj.species, in_obj.sp_taxid,idisolate)

    in_cursor.execute(query, val)


# dev:
#test_obj = ctables.info_16s()

#test_obj.isolate_name = "CC_PIBD_1_A2"
#test_obj.ab1_loc = "/home/xxxxxxxx"
#test_obj.sequence = "ATCG"
#test_obj.primer = "7f"

# goes to isolate table:
#test_obj.full_len_sequence = "ATCGAAAAAA"
#test_obj.lca_taxid = 816
#test_obj.sp_taxid = 816
#test_obj.species = "Clostridium"
