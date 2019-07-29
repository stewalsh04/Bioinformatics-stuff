# -*- coding: utf-8 -*-

import pandas as pd
import math

mydf=pd.read_csv('/home/steve/Myco/Myco_genome_features.csv', sep=',',header=None)

startofgenedf = pd.DataFrame()              #Set up dataframes
endofgenedf = pd.DataFrame()
centreofgenedf = pd.DataFrame()
lengthofgenes = pd.DataFrame()
intergene = pd.DataFrame()

mydf.columns = mydf.iloc[0]                 #Drop Top
mydf = mydf.drop([0,1,2], axis=0)
mydf2 = mydf

mydf.Start = pd.to_numeric(mydf.Start)      #Change to int64
mydf.End = pd.to_numeric(mydf.End)      

lengthofgenes['length']=mydf.End -mydf.Start   #Work out length of genes

###This section works out the starts and ends of genes

mydf['length']= lengthofgenes     #Add length to full df

startofgenedf['Chr'] = mydf.Accession    #This section creates start dataframe
startofgenedf['Start'] = mydf.Start - 50
startofgenedf['End'] = mydf.Start + ((mydf.length*0.1).astype(float).round())   #Takes first 10% of genes
startofgenedf['End'] = mydf.Start
startofgenedf['ID'] = mydf.Gene_ID + ".start"       #Adds start ID
startofgenedf['strand'] = mydf.Strand    #Adds strand

endofgenedf['Chr'] = mydf.Accession   ##Same as start for end
endofgenedf['Start'] = mydf.End - ((mydf.length*0.1).astype(float).round())
endofgenedf['End'] = mydf.End
endofgenedf['ID'] = mydf.Gene_ID + ".end"
endofgenedf['strand'] = mydf.Strand

centreofgenedf['Chr'] = mydf.Accession     #This section works out centre
centreofgenedf['Start'] = startofgenedf.End + 1    # uses outer regions to get centre
centreofgenedf['End'] = endofgenedf.Start -1
centreofgenedf['ID'] = mydf.Gene_ID + ".centre"
centreofgenedf['strand'] = mydf.Strand

frames = [centreofgenedf]   #Creates a list of what df's to merge
newdf = pd.concat(frames)  # Merges dfs

newdf = newdf.sort_values(by=['Start'])  #sorts new df by start


##This section works out intergene regions

nolast = mydf.drop(mydf.index[len(mydf)-1])
nolast = nolast.End
nolast = nolast.reset_index(drop=True)

mydf_skipfirstrow = mydf.iloc[1:]
mydf_skipfirstrow = mydf_skipfirstrow.Start
mydf_skipfirstrow = mydf_skipfirstrow.reset_index(drop=True)

intergene['chr']=mydf.Accession
intergene['Start']=nolast
intergene['End']=mydf_skipfirstrow
intergene['length']= intergene.End - intergene.Start
intergene = intergene.drop(intergene[intergene.length > 500].index)
intergene = intergene.drop(intergene[intergene.length < 5].index)

intergene = intergene.sort_values(by=['Start'])
intergene = intergene.reset_index(drop=True)

intergene['ID']= intergene.index + 1
intergene['ID']= intergene['ID'].apply(str)

intergene['ID'] = "INTER_" + intergene.ID
intergene = intergene.drop(intergene.index[len(intergene)-1,])
intergene = intergene.drop(intergene.index[len(intergene)-1,])
intergene = intergene.drop(intergene.index[len(intergene)-1,])
intergene = intergene.drop(intergene.index[len(intergene)-1,])

mydf2 = mydf2.drop("length", axis=1)
intergene = intergene.drop("length",axis=1)
print(startofgenedf)

startofgenedf.to_csv('/home/steve/Myco/starminus50.bed',index= False,sep="\t",float_format='%g')
