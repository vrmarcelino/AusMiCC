#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Helper functions related to MySQL

Created on Wed Jul 15 11:00:42 2020

@author: V.R.Marcelino
"""


import mysql.connector


### SQL settings
MYSQL_HOST='localhost'
MYSQL_USER='root'
MYSQL_PWD='4*4Genomes'
MYSQL_DB='ausmicc_001'


### Function to connect to the database
def db_connection():
    # from Sam's script
	#try to connect to the database otherwise die
	connection = None
	try:
		#SQLCONNECT: for adding sample and raw files
		connection = mysql.connector.connect(host=MYSQL_HOST, user=MYSQL_USER,
									 password=MYSQL_PWD, db=MYSQL_DB)
		return connection
    
	except:
		errmsg = "ERROR: Unable to connect to the database"
		print(errmsg)
		exit(-1)


### Function to get a tupple of the editable column names
### excludes the first (idxxx) and the last (timestamp) columns
def get_col_names(table_name, in_cursor):
    col_names_l=[]
    in_cursor.execute("SHOW COLUMNS FROM {};".format(table_name))

    for x in in_cursor:
        col_names_l.append(x[0])
    
    del col_names_l[0]# skip the first header (id, automatically generated)
    del col_names_l[-1] # skip the last header (timestamp, automatically generated)
    col_names = tuple(col_names_l)
    return (col_names)
    


### Function to find numeric identifiers for column names
### Requires a list of unique names, the table storing the unique ids and the cursor
### returns a dictionary with names and corresponding ids.
def get_ids(names, table, cursor):
    names2ids = {}   
    
    # guess table and idnames
    idname = "id" + table
    
    if table == 'pure_culture':
        col_name = 'ausmicc_name'
    else:
        col_name = table + "_name"

   
    for name in names:
        query="SELECT {} FROM {} WHERE {} = '{}'".format(idname, table, col_name, name)
    
        # try to fetch data, return an error if more than one id is found per sample
        cursor.execute(query)
        wanted_ids = cursor.fetchall()
        if len(wanted_ids) > 1:
            print ("more than one entry found for name %s" %(name))
            exit(-1)
        else:
            names2ids[name] = wanted_ids [0][0]
     
    return (names2ids)

 

