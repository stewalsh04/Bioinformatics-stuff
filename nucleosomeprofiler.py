# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 18:34:52 2019

@author: steve

This will work out the nucleosome profiles for a list of fragment frequencies such as that outputed from picard

"""
import numpy as np
import pandas as pd
import math
from pandas import DataFrame
import glob
import os

cwd = os.getcwd()

all_files = glob.glob(cwd + "/*.fraglengths") ##Ammend to reflect your file name, this is set for picard output

li2={}
totalnonenucleosome = {}
totalmononucleosome = {}
totaldinnuclesome = {}
totaltrinnuclesome = {}

bigtable = pd.DataFrame()

filenamesplit = ''
thisfilename = ''

for filename in all_files: ##This puts each file into a dict
    df = pd.read_csv(filename, index_col=None, sep='\t')
    filenamesplit = filename.split('.') ###This will need changing for your file name format!!!!!!!!!!!!!!!
    thisfilename = filenamesplit[0]
    thisfilename = thisfilename.split('-') ###This will need changing for your file name format!!!!!!!!!!!!!!!
    thisfilename = thisfilename[1]
    li2[thisfilename] = df

for i in li2: ##This extracts the nucleosome profiles and adds to a dict for each
    thislist = li2[i]
    nonenucleosome = li2[i]
    nonenucleosome = nonenucleosome[nonenucleosome['frag_length'] < 170]  
    nonenucleosome = nonenucleosome.sort_values(by=['frag_length'])
    totalnonenucleosome[i] = nonenucleosome['frequency'].sum()
    monnuclesome = thislist[thislist['frag_length'] > 170]  
    monnuclesome = monnuclesome[monnuclesome['frag_length'] < 300]
    monnuclesome = monnuclesome.sort_values(by=['frag_length'])
    totalmononucleosome[i] = monnuclesome['frequency'].sum()
    dinnuclesome = thislist[thislist['frag_length'] > 300]  
    dinnuclesome = dinnuclesome[dinnuclesome['frag_length'] < 500]
    dinnuclesome = dinnuclesome.sort_values(by=['frag_length'])
    totaldinnuclesome[i] = dinnuclesome['frequency'].sum()
    trinnuclesome = thislist[thislist['frag_length'] > 500]
    totaltrinnuclesome[i] = trinnuclesome['frequency'].sum()


none = pd.DataFrame.from_dict(totalnonenucleosome, orient='index') ##This converts the dicts to panda arrays
none.columns = ['none']
mono = pd.DataFrame.from_dict(totalmononucleosome, orient='index')
mono.columns = ['mono']
di = pd.DataFrame.from_dict(totaldinnuclesome, orient='index')
di.columns = ['di']
tri = pd.DataFrame.from_dict(totaltrinnuclesome, orient='index')
tri.columns = ['tri']

bigtable = none ##Next few line add individual panda arrays to a big array
bigtable['mono'] = mono.mono
bigtable['di'] = di.di
bigtable['tri'] = tri.tri
bigtable.loc['Total',:]= bigtable.sum(axis=0) ##Totals the columns


output = cwd + '/Nucleosome_profiles.csv'
bigtable.to_csv(output,index= True,sep="\t",float_format='%g')
