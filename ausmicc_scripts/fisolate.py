#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AusMiCC -  entry isolate info

Created on Wed Jul 15 13:10:33 2020

@author: V.R.Marcelino
"""


# local imports
from ausmicc_scripts import fconnector
from ausmicc_scripts import fparsedb


# debugging:
#sample_fp = args.input_fp
#df_fp = '/Users/vmar0011/Sync/02_Projects/0_AUS Microbiome/Example_entry_tables/Phase1_isolate_info.txt'

#in_df = pd.read_csv(df_fp, sep='\t')
#in_df = in_df.dropna(how='all')


def format_table(in_df):

    # convert yes/no to 1/0:
    in_df['spore_former'] = in_df['spore_former'].map(dict(yes=1, no=0))
    in_df['growth_plate'] = in_df['growth_plate'].map(dict(yes=1, no=0))
    in_df['growth_broth'] = in_df['growth_broth'].map(dict(yes=1, no=0))
    in_df['growth_aerobic'] = in_df['growth_aerobic'].map(dict(yes=1, no=0))
    in_df['growth_anaerobic'] = in_df['growth_anaerobic'].map(dict(yes=1, no=0))
    in_df['growth_microaerophilic'] = in_df['growth_microaerophilic'].map(dict(yes=1, no=0))


    # convert dates to a datetime obj, so it works with MySQL:
    in_df = fparsedb.conv2datetime(in_df,'isolation_date')
 
   
    # convert Nan to None, so it works with MySQL:
    in_df = fparsedb.convNaN2None(in_df)

    print ("Dataset formatted\n")
    return (in_df)




def add_plate(in_df, in_cursor):

    #### add plate information:
    
    plate_info = in_df[['plate_name','plate_freezer_location']].drop_duplicates()

    # check if every plate has a unique freezer location:
    if plate_info['plate_name'].is_unique == False:
        print ("""inconsistent plate and freezer location,
           the input table might have multiple freezer positions for the same plate,
           please correct""")
        exit(-1)
    
    for index, row in plate_info.iterrows():
        
        pl_name = row['plate_name']
        # check if plate exists in db first:
        query = "SELECT plate_name FROM plate WHERE plate_name = %s"
        val = (pl_name,)
        
        in_cursor.execute(query, val)
        existing_plate_name = in_cursor.fetchall()
        
        # if exists, skip it:
        if len(existing_plate_name) > 0:
            print ("Plate %s already in database, skipping." %(pl_name))
        
        else:
            query = "INSERT INTO plate (plate_name,plate_freezer_location) VALUES (%s, %s)"
            val = (row['plate_name'], row['plate_freezer_location'])
  
            in_cursor.execute(query, val)
            print ("Added plate %s to the database\n" %(pl_name))
        



def add_isolate(in_df, in_cursor):
    
    # get unique sample_names
    sample_names = in_df['sample_name'].unique()

    # use function to get a dict of sample_names:ids 
    sampleid_map = fconnector.get_ids(sample_names,'sample',in_cursor)

    # use function to get column names, drop teh ones we don't need in isolate table
    col_names = fconnector.get_col_names('isolate', in_cursor)

    # use function to get a dict of plate_names:ids
    plate_names = in_df['plate_name'].drop_duplicates()
    plateid_map = fconnector.get_ids(plate_names,'plate',in_cursor)


    ### add isolate information
    for index, row in in_df.iterrows():

        query = "INSERT INTO isolate {} VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(col_names).replace("'","`")
        val = (row['isolate_name'], row['isolate_barcode'],sampleid_map[row['sample_name']],
           plateid_map[row['plate_name']],None,None,row['isolation_date'], row['isolation_media'],
           row['colony_color'],row['colony_margin'],row['colony_form'],row['colony_elevation'],
           row['spore_former'],row['growth_plate'],row['growth_plate_image'],
           row['growth_broth'], row['growth_broth_image'], row['growth_aerobic'],
           row['growth_anaerobic'], row['growth_microaerophilic'], None, None, None)

        in_cursor.execute(query, val)

        ### Add isolate notes, when they exist
        if row['isolate_notes'] != None:
        
            # retrieve isolate id
            iso_name = [row['isolate_name']]
            isolateid_map = fconnector.get_ids(iso_name,'isolate',in_cursor)
        
            #add info to isolate_notes table
            query = "INSERT INTO isolate_notes (idisolate, note) VALUES (%s, %s)"
            val = (isolateid_map[row['isolate_name']], row['isolate_notes'])
        
            in_cursor.execute(query, val)
        
    print ("Isolate info added to database\n")




  
  
  
  
  
  
    