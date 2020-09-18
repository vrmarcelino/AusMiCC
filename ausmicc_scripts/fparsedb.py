#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Helper functions to parse the entry databases

Created on Wed Jul 15 11:11:09 2020

@author: V.R.Marcelino
"""

import pandas as pd



# convert dates to a datetime obj, so it works with MySQL:
def conv2datetime(in_df, col_name):
    # if 4 digits are provided for year, then %Y is required, otherwise we need a %y, so check which one it is:
    check_year_format = in_df[col_name][1].split("/")[-1]
    if len(check_year_format) == 2:
        in_df[col_name] = pd.to_datetime(in_df[col_name],format='%d/%m/%y')
    else:
        in_df[col_name] = pd.to_datetime(in_df[col_name], format='%d/%m/%Y')
    return (in_df)


# convert Nan to None, so it works with MySQL:
def convNaN2None(in_df):
    in_df = in_df.astype(object).where(pd.notnull(in_df), None)
    return (in_df)


# read files that will be saved as blob
# "rb" mode opens the file in binary format for reading, 
def read_bloob(filename):
    with open(filename, 'rb') as f:
        blob_file = f.read()
    return (blob_file)


