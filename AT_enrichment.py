# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

centre_df=pd.read_csv("/home/steve/Myco/centre.fasta", sep=',',header=None)
centre_df.columns = ['A']

start_df=pd.read_csv("/home/steve/Myco/start_end.fasta", sep=',',header=None)
start_df.columns = ['A']

ID = pd.DataFrame()
sequences = pd.DataFrame()


ID['ID'] = centre_df[centre_df['A'].str.contains('>')]
ID = ID[pd.notnull(ID['ID'])]
ID = ID.reset_index(drop=True)

sequences['sequence'] = centre_df[~centre_df['A'].str.contains('>')]
sequences = sequences[pd.notnull(sequences['sequence'])]
sequences = sequences.reset_index(drop=True)
Counts_df = pd.concat([ID, sequences], axis=1, join='inner')

Counts_df['sequence'] = Counts_df['sequence'].astype('|S')

Counts_df['AT'] = Counts_df['sequence'].str.count('A') + Counts_df['sequence'].str.count('T')
Counts_df['GC'] = Counts_df['sequence'].str.count('G') + Counts_df['sequence'].str.count('C')
Counts_df['AT%'] = (Counts_df['AT'] / (Counts_df['AT'] + Counts_df['GC'])) * 100
Counts_df['GC%'] = (Counts_df['GC'] / (Counts_df['GC'] + Counts_df['AT'])) * 100


ID2 = pd.DataFrame()
sequences2 = pd.DataFrame()


ID2['ID'] = start_df[start_df['A'].str.contains('>')]
ID2 = ID2[pd.notnull(ID2['ID'])]
ID2 = ID2.reset_index(drop=True)

sequences2['sequence'] = start_df[~start_df['A'].str.contains('>')]
sequences2 = sequences2[pd.notnull(sequences2['sequence'])]
sequences2 = sequences2.reset_index(drop=True)
Counts_df2 = pd.concat([ID2, sequences2], axis=1, join='inner')

Counts_df2['sequence'] = Counts_df2['sequence'].astype('|S')

Counts_df2['AT'] = Counts_df2['sequence'].str.count('A') + Counts_df2['sequence'].str.count('T')
Counts_df2['GC'] = Counts_df2['sequence'].str.count('G') + Counts_df2['sequence'].str.count('C')
Counts_df2['AT%'] = (Counts_df2['AT'] / (Counts_df2['AT'] + Counts_df2['GC'])) * 100
Counts_df2['GC%'] = (Counts_df2['GC'] / (Counts_df2['GC'] + Counts_df2['AT'])) * 100

#Total_AT_centre = (Counts_df['AT'].sum()) / ((Counts_df['AT'].sum()) + (Counts_df['GC'].sum())) #* 100
Total_AT_centre = Counts_df['AT'].sum() 
Totalbases_centre = Counts_df['AT'].sum() + Counts_df['GC'].sum()
percentage_AT_centre = np.true_divide (Total_AT_centre, Totalbases_centre) * 100
Total_AT_start = Counts_df2['AT'].sum()
Totalbases_start = Counts_df2['AT'].sum() + Counts_df2['GC'].sum()
percentage_AT_start = np.true_divide (Total_AT_start, Totalbases_start) * 100

print(Counts_df)
print(Counts_df2)
print(percentage_AT_centre)
print(percentage_AT_start)
