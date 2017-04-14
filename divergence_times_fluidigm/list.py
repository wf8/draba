#! /usr/bin/python

input_dir = "data/org_fluidigm_data_fasta/"

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from os import listdir
from os.path import isfile, join

files = [ f for f in listdir(input_dir) if isfile(join(input_dir,f)) ]

# samples to omit
skip = []

# get list of all samples
taxa = []
for f in files:
    handle = open(input_dir + f, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        description = str(record.description)
        if description not in taxa and description not in skip:
            taxa.append(description)

print("Total num samples = " + str(len(taxa)))
for x in sorted(taxa):
    print(x)

