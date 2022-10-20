#!/usr/bin/env python
## Copyright (c) 2014 Marnix H. Medema
## Max Planck Institute for Marine Microbiology
## Microbial Genomics & Bioinformatics Research Group

###Script that extract trans-KS domains out of protein sequences

#Libraries to import
import sys
import os
import string
import subprocess

try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

import warnings
with warnings.catch_warnings(): # Don't display the SearchIO experimental warning.
    warnings.simplefilter("ignore")
    from Bio import SearchIO


#Global settings
#global cpu
##Number of CPU cores used for multiprocessing
#cpu = 4

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
    except OSError, e:
        print "%r %r returned %r" % (commands, input[:40], e)
        raise

def run_hmmsearch(query_hmmfile, fasta_file):
    "Run hmmsearch"
    command = ["hmmsearch", "--cpu", "3","--domT", "50", "--domtblout",   			"temp.out", query_hmmfile, fasta_file]
    try:
        out, err, retcode = execute(command)
	out = open("temp.out","r").read()
    except OSError:
        return []
    res_stream = StringIO(out)
    results = list(SearchIO.parse(res_stream, 'hmmsearch3-domtab'))
    return results

def read_fasta_file(fastafile):
    "Get sequences from fasta file and store in dictionary"
    ###My idea:
    infile = open(fastafile).read()
    entries = infile.split(">")[1:]
    #print entries
    fastadict = {} #Accession as key, sequence as value
    for entry in entries:
        accession = entry.partition("\n")[0].replace("\r","")
        #print accession
        sequence = entry.partition("\n")[2].partition("\n")[0]
        #print sequence
        fastadict[accession] = sequence  
    return fastadict

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
                        print domain_name, "= cis-AT; AT domain score:", at_hit.bitscore, "From position", at_hit.hit_start, "to", at_hit.hit_end
                        break
            #If not, write sequence to output file
            if transAT and len(domain_sequence) > 380: ##Don't get truncated KS domains
                out_file.write(">%s\n%s\n" % (domain_name, domain_sequence))
            domnr += 1
        domnr = 1
    out_file.close()


if __name__ == "__main__":
    "Run this file from here as a script"

    #Read command-line parameters
    ks_hmmfile = "PKS_KS.hmm"
    at_hmmfile = "PKS_AT.hmm"
    fastafile = "knownclusterprots.fasta"
    outfile = "transKS_domains.fasta"
    
    #Read FASTA file
    fastadict = read_fasta_file(fastafile)

    #Run HMMer
    ks_hmmer_results = run_HMMer(ks_hmmfile, fastafile)
    at_hmmer_results = run_HMMer(at_hmmfile, fastafile)

    #Write FASTA file with domains
    write_fasta(fastadict, outfile, ks_hmmer_results, at_hmmer_results)

