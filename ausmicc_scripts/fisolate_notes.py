#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AusMiCC -  entry isolate_notes

Created on Thu Jul 16 10:39:18 2020

@author: V.R.Marcelino
"""

from ausmicc_scripts import fconnector


# function will take filepath and in_isolate_name

def add_isolate_notes(in_df, in_cursor,in_isolate_name):

    # idenitfy the idisolate
    iso_name = [in_isolate_name]

    sampleid_map = fconnector.get_ids(iso_name,'isolate',in_cursor)

    query = "INSERT INTO isolate_notes (idisolate,note) VALUES (%s,%s)"
    val = (sampleid_map[iso_name[0]], in_df)

    in_cursor.execute(query, val)

    print ("\nIsolate_notes added to database")
    
    
