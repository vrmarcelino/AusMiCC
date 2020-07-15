#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to enter data into AusMiCC. Staring with info from phase 1.

Using a local copy of the database

Created on Thu Jul  2 08:34:16 2020

@author: V.R.Marcelino
"""

#imports
import pandas as pd
from argparse import ArgumentParser


# local imports
from ausmicc_scripts import fconnector
from ausmicc_scripts import fsample
from ausmicc_scripts import fisolate

### input files
parser = ArgumentParser()
parser.add_argument('-i', '--input_fp', help="""Path to the tab_delimited input file""", required=True)
parser.add_argumment ('-m', '--mode', help = """ Which data type are you entering? 
                      Options: sample, isolate, isolate_notes, etc""", required=True)


args = parser.parse_args()
df_fp = args.input_fp
mode = args.mode

#df_fp = '/Users/vmar0011/Sync/02_Projects/0_AUS Microbiome/Example_entry_tables/Phase1_isolate_info.txt'
#mode = 'isolate'


### connect to the database:
aus_db_conn = fconnector.db_connection()
cursor = aus_db_conn.cursor()



### Read the sample table
df = pd.read_csv(df_fp, sep='\t')

### delete rows that contain only NaNs
df = df.dropna(how='all')


### process sample info, use the sample function to format the table and add to the database:
if mode == 'sample':
    fsample.add_sample(df, cursor)
    
    aus_db_conn.commit()
    
    print ("\nDone! Sample data entered in the AusMiCC database. \n")
    print ("You can proceed to the next stage and enter isolate information. \n")



if mode == 'isolate':
    df = fisolate.format_table(df)  

    fisolate.add_plate(df, cursor)
  
    fisolate.add_isolate(df,cursor)

    aus_db_conn.commit()
    
    


       
#aus_db_conn.commit()
aus_db_conn.close()



# Notes:
#cursor.execute("SELECT * FROM isolate;")
#cursor.execute("DELETE FROM plate;")
#for x in cursor:
#  print(x)


#print(cursor.rowcount, "records")
#cursor.reset()




