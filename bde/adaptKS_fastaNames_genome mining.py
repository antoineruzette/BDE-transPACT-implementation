#!/usr/bin/env python
# coding: utf-8

# In[173]:


#the goal of that script is to rename the KS names in the fasta to match with the KS names in the clades.tsv file


# In[188]:


#load libraries
import pandas as pd
import numpy as np
import regex as re
import os
import sys

#change your working directory
os.chdir("C:\\Users\\antoi\\transPACT\\data\\genome_mining")
current_dir = os.getcwd()
os.chdir(current_dir)


# In[2]:


#not used anymore in genome mining approach
def KSnames_full2unique(df): 
    #improvement to code: check if the name is unique -> if not, incorporate the second name of the KS sequence
    specific_names = []
    unique_names = []
    
    for index, row in df.iterrows(): 
        #row[0] corresponds to the names without the df index
        line_name = re.split('\|', row[0]) #need to escape the separator with a backslash
        specific_names.append(line_name)
        
    remove_df = pd.DataFrame(specific_names)
    
    for lines in specific_names: 
        tmp = str(lines[0]) + '_' + str(lines[1])#+ '_' + str(lines[2])
        unique_names.append(tmp)
    unique_names = pd.DataFrame(unique_names)    
        
    return unique_names


# In[184]:


def extract_KS_info(file): 
    
    sequences = []
    names = []
    refseq = []
    pathwayID = []
    moduleID = []
    
    numbers = []
    pattern = re.compile(r'_(KS[0-9]{1,2})')
    
    df_trimmed = pd.DataFrame()
    
    with open(file) as infile:
        for lines in infile: 
            sequence_lines = re.split('^(>)|\n', lines)
            sequence_lines[-1] = sequence_lines[-1].replace('\n', '') #remove '\n' from the previous split
            #print(sequence_lines)
            sequences.append(sequence_lines[0])
            names.append(sequence_lines[2])
            
        #use list comprehension to remove empty elements    
        full_sequences = [element for element in sequences if element.strip()]
        #print(len(full_sequences))#correct -> 150
        full_names = [element for element in names if element.strip()]
        #print(len(full_names))#correct -> 150
        
        for name in full_names: 
            lines = re.split('\|', name)
            refseq.append(lines[0])
            pathwayID.append(lines[1])
            moduleID.append(lines[2])           
        #print(len(refseq))
        
        
        #extract only the KS numbers
        full_names = pd.Series(full_names) #need to convert to Series to use findall() function
        for name in full_names: 
            for match in re.findall(pattern, name): 
                numbers.append(match)
                
        #full_names = pd.DataFrame(full_names)
        #full_unique_names = pd.DataFrame(KSnames_full2unique(full_names))
        full_numbers = pd.DataFrame(numbers)
        full_sequences = pd.DataFrame(full_sequences)
        full_refseq = pd.DataFrame(refseq)
        full_pathwayID = pd.DataFrame(pathwayID)
        full_moduleID = pd.DataFrame(moduleID)
        
        #saving the extracted results in a main df 
        df = full_refseq.join(full_pathwayID, how = 'left', lsuffix = 'refseqID', rsuffix = 'pathwayID')
        df = df.join(full_moduleID, how = 'left', lsuffix='moduleID', rsuffix='moduleID')
        df = df.join(full_numbers, how = 'left', lsuffix='left', rsuffix='right')
        df = df.join(full_sequences, how = 'left', lsuffix='left', rsuffix='right')
        mapping = {df.columns[0]:'refSeqID', df.columns[1]: 'pathwayID', 
                   df.columns[2]:'moduleID', df.columns[3]: 'KS', df.columns[4]:'sequences'}
        df = df.rename(columns = mapping)

        #trimming df based on the number of KS domains making up a pathway
        #idea: reducing computing time by removing pathways adding up little relevant info to the analysis
        df_trimmed = df.groupby(by = ['refSeqID', 'pathwayID']).filter(lambda x: len(x) > 10 and len(x) < 16)

    return df, df_trimmed


# In[153]:


#df_pks = df_trimmed.groupby(['refSeqID', 'pathwayID'])
#df_pks = df_trimmed.groupby('refSeqID')
#test = df_pks.get_group('NZ_CP047644')
#test
#df_test.get_group('c1454370-1542202')


# In[183]:


def order_rename(df): 
    
    #idea
    df_global = pd.DataFrame()
    genome_list = pd.unique(df['refSeqID'])
    
    for genome in genome_list: 
        df_genome = df.loc[df['refSeqID'] == genome]
        pathway_list = pd.unique(df_genome['pathwayID'])
        
        for pathway in pathway_list:
            df_pathway = df_genome.loc[df_genome['pathwayID'] == pathway]
            df_pathway = df_pathway.sort_values(by = 'moduleID')
            
            for ks in range(0, len(df_pathway)): 
                ks_string = str(ks + 1)
                string = 'KS' + ks_string
                df_pathway.iloc[[ks],[3]] = string
            df_global = df_global.append(df_pathway)
            
    df_global = df_global.sort_values(by = ['refSeqID', 'pathwayID'], ascending = [True, True])
    return df_global


# In[208]:


annotation_file = 'tKS_domains.fasta'

df, df_trimmed = extract_KS_info(annotation_file)
df_trimmed_ordered = order_rename(df_trimmed)


# In[215]:


#check number of KSs and PKSs under study
print(len(df))
print(len(pd.unique(df['pathwayID'])))
df_test = df.groupby(by = ['refSeqID', 'pathwayID']).filter(lambda x: len(x) > 2)
print(len(df_test))
print(len(pd.unique(df_test['pathwayID'])))


# In[191]:


with open('full_trimmed_renamed.fasta', 'w') as outfile: 
    for row in df_trimmed_ordered.values: #for each row
        outfile.write('>' + str(row[0]) + '_' + str(row[1]) + '|' + str(row[3])
                      + '\n' + str(row[4]) + '\n') #print each column of each row
    outfile.close()

