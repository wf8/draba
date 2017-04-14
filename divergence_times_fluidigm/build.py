#! /usr/bin/python

#
# Script to construct Draba fluidigm sequence data matrices.
#

input_dir = "data/org_fluidigm_data_fasta/"
output_dir = "data/final_fluidigm_data_fasta/"

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from os import listdir
from os.path import isfile, join

files = [ f for f in listdir(input_dir) if isfile(join(input_dir,f)) ]

# samples to omit
skip = [
"Arabidopsis1_control",
"Arabidopsis2_control",
"Arabidopsis3_control",
"Arabidopsis4_control",
"Arabidopsis5_control",
"Arabidopsis6_control",
"Arabidopsis8_control",
"Dra_1002_D_ramosissima",
"Dra_1002_D_ramosissima_2",
"Draba_abajoensis_UC9_cld2",
"Draba_aizoides_UC25_cld1",
"Draba_daviesiae_UC13_cld2",
"Draba_novolympica_IEJT235_cld2",
"Draba_oligo_1310_Pop2016",
"Draba_oligo_1311_Pop2016",
"Draba_oligo_1312_Pop2016",
"Draba_oligo_1313_Pop2016",
"Draba_oligo_1315_Pop2016",
"Draba_oligo_1316_Pop2016",
"Draba_oligo_1317_Pop2016",
"Draba_oligo_1318_Pop2016",
"Draba_oligo_1319_Pop2016",
"Draba_oligo_1320_Pop2016",
"Draba_oligo_1321_Pop2016",
"Draba_oligo_1322_Pop2016",
"Draba_oligo_1323_Pop2016",
"Draba_oligo_1324_Pop2016",
"Draba_oligo_1325_Pop2016",
"Draba_oligo_1326_Pop2016",
"Draba_oligo_1327_Pop2016",
"Draba_oligo_1328_Pop2016",
"Draba_oligo_1329_Pop2016",
"Draba_oligo_1330_Pop2016",
"Draba_oligo_1331_Pop2016",
"Draba_oligo_1332_Pop2016",
"Draba_oligo_1333_Pop2016",
"Draba_oligo_1624_Pop2027",
"Draba_oligo_1625_Pop2027",
"Draba_oligo_1626_Pop2027",
"Draba_oligo_1627_Pop2027",
"Draba_oligo_1628_Pop2027",
"Draba_oligo_1629_Pop2027",
"Draba_oligo_1630_Pop2027",
"Draba_oligo_1631_Pop2027",
"Draba_oligo_1632_Pop2027",
"Draba_oligo_1633_Pop2027",
"Draba_oligo_1634_Pop2027",
"Draba_oligo_1635_Pop2027",
"Draba_oligo_1636_Pop2027",
"Draba_oligo_1637_Pop2027",
"Draba_oligo_1638_Pop2027",
"Draba_oligo_1640_Pop2027",
"Draba_oligo_1641_Pop2027",
"Draba_oligo_1642_Pop2027",
"Draba_oligo_1643_Pop2027",
"Draba_oligo_1644_Pop2027",
"Draba_oligo_1645_Pop2027",
"Draba_oligo_1646_Pop2027",
"Draba_oligo_1647_Pop2027",
"Draba_oligo_1648_Pop2027",
"Draba_oligo_1650_Pop2027",
"Draba_oligo_1651_Pop2027",
"Draba_oligo_1652_Pop2027",
"Draba_oligosperma_IEJT274_cld3",
"EKY_01_D_ramosissima",
"EKY_03_D_ramosissima",
"EKY_04_D_ramosissima",
"WVA_02_D_ramosissima",
"WVA_05_D_ramosissima",
"WVA_08_D_ramosissima",
"WVA_09_D_ramosissima",
"WVA_11_D_ramosissima",
"WVA_12_D_ramosissima",
"WVA_13_D_ramosissima",
"WVA_20_D_ramosissima"
        ]

# get list of all samples we are keeping
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

# now add missing sequences
for f in files:
    found_taxa = []
    final_records = []
    handle = open(input_dir + f, "rU")
    seq_length = 0
    for record in SeqIO.parse(handle, "fasta"):
        seq_length = len(str(record.seq))
        binomial = str(record.description)
        if binomial not in found_taxa and binomial not in skip:
            found_taxa.append(binomial)
        if binomial not in skip:
            final_records.append(record)
    # make missing sequences
    for taxon in taxa:
        if taxon not in found_taxa:
            final_records.append(SeqRecord(Seq("?" * seq_length), id=taxon, description=""))
    output_handle = open(output_dir + f , "w")
    SeqIO.write(final_records, output_handle, "fasta")
