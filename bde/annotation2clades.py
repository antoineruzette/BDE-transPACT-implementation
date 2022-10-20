#!/usr/bin/env python
# coding: utf-8

# In[22]:


#load libraries
import pandas as pd
import numpy as np
import regex as re
import os
import sys

#change your working directory
os.chdir("C:\\Users\\antoi\\transPACT\\data\\dendrogram_647KS")
current_dir = os.getcwd()
os.chdir(current_dir)


# In[23]:


def extract_KS_numbers(file): 
    "Extract the KS numbers in the order of the file and store them in a list"
    
    KS_numbers = []
    
    with open(file) as infile:
        next(infile) #skip header
        for lines in infile: 
            pattern = re.compile(r'(KS\d+)')
            for match in re.findall(pattern, lines):
                #print(match)
                KS_numbers.append(match)
        KS_numbers = pd.DataFrame(KS_numbers)
    return KS_numbers


# In[24]:


def extract_KS_desc(file): 
    
    KS_desc = []
    
    with open(file) as infile: 
        next(infile)#skip header
        for lines in infile: 
            lines_split = re.split('\t', lines)
            lines_split[-1] = lines_split[-1].replace('\n', '')
            KS_desc.append(lines_split[-1])
            
        KS_desc = pd.DataFrame(KS_desc)
    return KS_desc


# In[25]:


def KSnames_full2unique(df): 
    #improvement: check if the name is unique -> if not, incorporate the second name of the KS sequence
    specific_names = []
    
    for index, row in df.iterrows(): 
        line_name = re.split('(_)', row[0])
        #line_name = lines_name.replace('\n', '')
        specific_names.append(line_name)
        #print(row[0])#corresponds to the names without the df index
    specific_names = pd.DataFrame(specific_names)
        
    KS_unique_names = (pd.DataFrame(specific_names[0])).astype(str)
    return KS_unique_names


# In[26]:


def extract_KS_info(file): 
    
    entries = []
    
    with open(file, 'r') as infile:
        next(infile) #skip header
        for lines in infile: 
            lines_split = re.split('(_KS\d+[_])|\t', lines) #split by '_KS[digit]_' and 'tab \t'
            lines_split[-1] = lines_split[-1].replace('\n', '') #remove '\n' from the previous split
            
            #lines_split[1] = lines_split[1].replace('_', '') #remove '_' from the previous split 
            entries.append(lines_split)            
            #print(lines_split)
        
        entries = pd.DataFrame(entries)
        #KS_full_names = (pd.DataFrame(entries[0] + entries).astype(str)
        KS_full_names = (pd.DataFrame(entries[0] + entries[1] + entries[2])).astype(str)
        KS_clade_numbers = (pd.DataFrame(entries[4])).astype(str)
        KS_numbers = (extract_KS_numbers(file)).astype(str)
        KS_desc = (extract_KS_desc(file)).astype(str)
        
        KS_unique_names = pd.DataFrame(KSnames_full2unique(KS_full_names))
            
    return KS_unique_names, KS_clade_numbers, KS_numbers, KS_desc


# In[27]:


#read files

annotation_file = 'annotation_ar_20220331.tsv'
#needed to remove the _KS0_ from calyculin_CalB_1_KS3_KS0_b_L_OH at ln.185 from the annotation file (annotation_mc_20180727.tsv)
#otherwise, the number of rows do not match between the different extractions

#call functions
names, clades, numbers, desc = extract_KS_info(annotation_file)
df = names.join(clades) #merge into a single df
df = df.join(numbers, how='left', lsuffix='KS_name')
df = df.join(desc, how = 'left', rsuffix = 'KS_desc')
df = df.rename(columns = {'0KS_name':'KS_full_names', 4: 'KS_clade', '0': 'KS_number'})


# In[28]:


#write the annotation.clades.tsv file: PKS|KS \t Clade_x 

with open('full_annotation.clades.tsv', 'w') as outfile: 
    for row in df.values: #for each row
        outfile.write(row[0] + '|' + row[2] + '\t' + row[1] + '\n') #print each column of each row
    outfile.close()


# In[29]:


#write the renamed annotation file with PKS|KS \t Clade_x \t Clade_desc

with open('full_annotation.tsv', 'w') as outfile: 
    for row in df.values: #for each row
        outfile.write(row[0] + '|' + row[2] + '\t' + row[1] + '\t' + row[3] + '\n') 
    outfile.close()

