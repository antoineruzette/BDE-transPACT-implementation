#!/usr/bin/env python
# coding: utf-8

# In[176]:


#the goal of that script is to rename the KS names in the fasta to match with the KS names in the clades.tsv file


# In[239]:


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


# In[240]:


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


# In[241]:


def extract_KS_numbers(df_names): 
    "Extract the KS numbers in the order of the file and store them in a list"
    
    #df_names = pd.DataFrame(df_names)
    KS_numbers = []
    pattern = re.compile(r'(_KS\d+_)')
    
    for names in df_names.iterrows(): 
        print(names)
        #for match in re.findall(pattern, names[0]): 
            #KS_numbers.append(match)
    KS_numbers = pd.DataFrame(KS_numbers)
    return KS_numbers


# In[242]:


def KSnames_full2unique(df): 
    #improvement to code: check if the name is unique -> if not, incorporate the second name of the KS sequence
    specific_names = []
    
    for index, row in df.iterrows(): 
        line_name = re.split('(_)', row[0])
        #line_name = lines_name.replace('\n', '')
        specific_names.append(line_name)
        #print(row[0])#corresponds to the names without the df index
    specific_names = pd.DataFrame(specific_names)
        
    KS_unique_names = (pd.DataFrame(specific_names[0])).astype(str)
    return KS_unique_names


# In[243]:


def extract_KS_info(file): 
    
    sequences = []
    names = []
    numbers = []
    pattern = re.compile(r'_(KS\d+)_')
    
    with open(file) as infile:
        for lines in infile: 
            sequence_lines = re.split('^(>)|\n', lines)
            sequence_lines[-1] = sequence_lines[-1].replace('\n', '') #remove '\n' from the previous split
            sequences.append(sequence_lines[0])
            names.append(sequence_lines[2])
            
        #use list comprehension to remove empty elements    
        full_sequences = [element for element in sequences if element.strip()]
        full_names = [element for element in names if element.strip()]
        
        #extract only the KS numbers
        full_names = pd.Series(full_names)#need to convert to Series to use findall() function
        for name in full_names: 
            for match in re.findall(pattern, name): 
                numbers.append(match)
        
        full_names = pd.DataFrame(full_names)
        full_unique_names = pd.DataFrame(KSnames_full2unique(full_names))
        full_numbers = pd.DataFrame(numbers)
        full_sequences = pd.DataFrame(full_sequences)     
        #KS_numbers = extract_KS_numbers(full_names)
    return full_sequences, full_unique_names, full_numbers


# In[244]:


annotation_file = '647KS_sequences.fasta'

sequences, names, numbers = extract_KS_info(annotation_file)
if len(sequences) == len(names) == len(numbers): #length of numbers is 648
    print(f'Each df length match to', len(sequences))
else: print(f'len_KS_numbers', len(numbers))

df = names.join(sequences, how = 'left', lsuffix = 'names', rsuffix = 'sequences')
df = df.join(numbers)
df

#problems with the following KS sequences: do not possess a KS number
#chlorotonil_CtoD_4_aMeshDB -> added KS8 in its name
#rhizoxins_RhiD_1_K9_aMeb_D_OH -> corrected the K9 into KS9
#2qo3_chainA_EryKS3_KS1_OUTGROUP -> added KS1 to the two outgroups
#bryostatin_BryX_3_KS16_added -> added _added at the end (so that the pattern can be matched)


# In[245]:


with open('647KS_renamed.fasta', 'w') as outfile: 
    for row in df.values: #for each row
        outfile.write('>' + str(row[0]) + '|' + str(row[2]) + '\n' + str(row[1]) + '\n') #print each column of each row
    outfile.close()

