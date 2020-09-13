#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AusMiCC -  entry sample info

Created on Wed Jul 15 11:36:35 2020

@author: V.R.Marcelino
"""

# local imports
from ausmicc_scripts import fconnector
from ausmicc_scripts import fparsedb



def add_sample(in_df, in_cursor):
    
    # convert dates to a datetime obj, so it works with MySQL:
    in_df = fparsedb.conv2datetime(in_df,'sample_date')
    in_df = fparsedb.conv2datetime(in_df,'pat_birth_date')

    # convert Nan to None, so it works with MySQL:
    in_df = fparsedb.convNaN2None(in_df)

    #### get a list of column names:
    col_names = fconnector.get_col_names('sample',in_cursor)


    #### add to the mysql database, line by line:
    for index, row in in_df.iterrows():
        query = "INSERT INTO sample {} VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)".format(col_names).replace("'","`")
        val = (row['sample_name'],row['project_name'],row['researcher_name'],row['body_site'],
               row['sample_type'],row['sample_date'],row['pat_age_years'],
               row['pat_birth_date'],row['city'],row['state'],row['country'],row['ethnicity'])
    
        in_cursor.execute(query, val)
        
    print ("query complete")
    
    
    
    

