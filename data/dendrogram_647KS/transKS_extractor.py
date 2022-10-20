#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Libraries to import
import sys
import os
import string
import subprocess
import warnings
import regex as re
import pandas as pd
import numpy as np


try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


# In[94]:


#Global settings and working directory
global cpu
##Number of CPU cores used for multiprocessing
cpu = 4

os.chdir("C:\\Users\\antoi\\transPACT")
os.getcwd()


# In[95]:


def execute(commands, input=None):
    "Execute commands in a system-independent manner"

    if input is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    proc = subprocess.Popen(commands, stdin=stdin_redir,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    try:
        out, err = proc.communicate(input=input)
        retcode = proc.returncode
        return out, err, retcode
    except (OSError, e):
        print("%r %r returned %r" % (commands, input[:40], e))
        raise


# In[96]:


def run_hmmsearch(query_hmmfile, fasta_file):
    "Run hmmsearch"
    command = ["hmmsearch", "--domT", "50", "--domtblout", "temp.out", query_hmmfile, fasta_file]
    try:
        out, err, retcode = execute(command)
        out = open("temp.out","r").read()
    except OSError:
        return []
    res_stream = StringIO(out)
    results = list(SearchIO.parse(res_stream, 'hmmsearch3-domtab'))
    return results


# In[100]:


def read_fasta_file(fastafile):
    "Get sequences from fasta file and store in dictionary"
    #adapted by Antoine Ruzette
    entries = []
    names = []
    sequences = []
    with open(fastafile) as infile: 
        for lines in infile: 
            lines_split = re.split('(^>)', lines)
            entries.append(lines_split)
        entries = pd.DataFrame(entries)
        entries = entries.rename(columns = { 2:'names', 1: 'rest', 0: 'sequences'})

        #split names and sequences in two 
        names = pd.DataFrame(entries['names'])
        names = pd.DataFrame((names.replace('', np.nan)).dropna())
        names = names.reset_index(drop = True)
        sequences = pd.DataFrame(entries['sequences'])
        sequences = pd.DataFrame((sequences.replace('', np.nan)).dropna())
        sequences = sequences.reset_index(drop = True)
        
        #bring names and sequences back together
        entries = pd.DataFrame()
        entries['names'] = names
        entries['sequences'] = sequences
        entries = entries.replace('\n', '')
        
        #sent it to dictionary
        
        entries = dict(entries.values)
        entries = {key.strip(): item.strip() for key, item in entries.items()}#remove special caracter \n
            
    return entries, names, sequences


# In[101]:


#HMMer is a sequence analysis package
def run_HMMer(hmmfile, fastafile):
    "Run HMMer to find start and end sites of KS/AT domains"
    runresults = run_hmmsearch(hmmfile, fastafile)
    results_by_id = {}
    for runresult in runresults:
        #Store results in dictionary by NRPS accession
        for hsp in runresult.hsps:
            if not results_by_id.has_key(hsp.hit_id):
                results_by_id[hsp.hit_id] = [hsp]
            else:
                results_by_id[hsp.hit_id].append(hsp)
    return results_by_id


# In[102]:


def write_fasta(fastadict, outfile, ks_hmmer_results, at_hmmer_results):
    #For each HMM, print a FASTA file with all domain hits
    out_file = open(outfile, "w")
    domnr = 1
    for cds in ks_hmmer_results.keys():
        #Get sequence from fastadicts
        cds_sequence = fastadict[cds]
        #For each hit, write out domain name + sequence in FASTA format
        for hit in ks_hmmer_results[cds]:
            domain_name = "%s_KS%s" % (cds, domnr)
            domain_sequence = cds_sequence[hit.hit_start:hit.hit_end]
            #Check if there is no AT domain in the 400 AA after this KS domain
            transAT = True
            if at_hmmer_results.has_key(cds):
                at_hits = at_hmmer_results[cds]
                for at_hit in at_hits:
                    if at_hit.hit_start > hit.hit_end and int(at_hit.hit_start) - int(hit.hit_end) < 400:
                        transAT = False
                        print(domain_name, "= cis-AT; AT domain score:", at_hit.bitscore, "From position", at_hit.hit_start, "to", at_hit.hit_end)
                        break
            #If not, write sequence to output file
            if transAT and len(domain_sequence) > 380: ##Don't get truncated KS domains
                out_file.write(">%s\n%s\n" % (domain_name, domain_sequence))
            domnr += 1
        domnr = 1
    out_file.close()


# In[103]:


#Read command-line parameters
#Replace the file names by yours
ks_hmmfile = "PKS_KS.hmm"
at_hmmfile = "PKS_AT.hmm"
#fastafile1 = "subset_aligned.afa"#aligned fasta file
fastafile2 = "knownclusterprots.fasta"
outfile = "transKS_domains2.fasta"

#reading the fasta file works => sent to a dictionary with names as key and sequences as content
entries, names, sequences = read_fasta_file(fastafile2)
#.hmm files contain the profile hmm architecture. How was it generated? => use hmmbuild function on a multiple sequence alignment


# In[106]:


#Run HMMer
ks_hmmer_results = run_HMMer(ks_hmmfile, fastafile2)
at_hmmer_results = run_HMMer(at_hmmfile, fastafile2)
   
#Write FASTA file with domains
write_fasta(entries, outfile, ks_hmmer_results, at_hmmer_results)

