#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to enter data into AusMiCC. Staring with info from phase 1.

Using a local copy of the database

Created on Thu Jul  2 08:34:16 2020

@author: V.R.Marcelino
"""

#imports
import mysql.connector # not used in this script, but needed to avoid a segmentation fault
import pandas as pd
from argparse import ArgumentParser
import sys


# local imports
from ausmicc_scripts import fconnector
from ausmicc_scripts import fparsedb
from ausmicc_scripts import fsample
from ausmicc_scripts import fisolate
from ausmicc_scripts import fisolate_notes

### input files
parser = ArgumentParser()

parser.add_argument('-i', '--input_fp', help="""Path to the tab_delimited input file 
                    (other file formats are acceptable for isolate_notes only)""", required=True)
parser.add_argument('-m', '--mode', help = """ Which data type are you entering? 
                      Options: sample, isolate, isolate_notes, etc""", required=True)
parser.add_argument('-in', '--isolate_name', help = """Name of the isolate 
                      to which these notes refers to (only required when adding isolate_notes)""", required=False)

args = parser.parse_args()
df_fp = args.input_fp
mode = args.mode
in_isolate_name = args.isolate_name

#df_fp = '/home/vmar0011/Phase1_isolate_info_ab1_test.txt'
#mode = 'isolate'
#in_isolate_name = "AMR001_A2"

### connect to the database:
aus_db_conn = fconnector.db_connection()
cursor = aus_db_conn.cursor()


### If this is a note (bloob), must be treated differently:
if mode == "isolate_notes":

    # open isolate_notes as binary
    bloob_file = fparsedb.read_bloob(df_fp)
    
    # save note
    fisolate_notes.add_isolate_notes(bloob_file,cursor,in_isolate_name)
    
    aus_db_conn.commit()
    sys.exit("Done!\n")
    


# Read the sample table
df = pd.read_csv(df_fp, sep='\t')
    
# delete rows that contain only NaNs
df = df.dropna(how='all')

### process sample info, use the sample function to format the table and add to the database:
if mode == 'sample':
    fsample.add_sample(df, cursor)
    
    aus_db_conn.commit()
    
    print ("\nDone! Sample data entered in the AusMiCC database. \n")
    print ("You can proceed to the next stage and enter isolate information. \n")



if mode == 'isolate':
    df = fisolate.format_table(df)  

    fisolate.add_plate(df,cursor)
  
    fisolate.add_isolate(df,cursor)

    aus_db_conn.commit() 
    print ("\nDone!\n")

#aus_db_conn.commit()
aus_db_conn.close()



# Notes:
#cursor.execute("SELECT * FROM isolate;")
#cursor.execute("DELETE FROM plate;")
#for x in cursor:
#  print(x)


#print(cursor.rowcount, "records")
#cursor.reset()




