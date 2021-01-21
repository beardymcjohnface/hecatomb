"""

Just download the databases. You should only need to do this once (unless we update the databases!)

Rob Edwards, Feb 2020
Updated: Michael Roach, Jan 2021
"""

import os
import sys


# load config file
if not config:
    sys.stderr.write("FATAL: Please define a config file using the --configfile command line option.\n")
    sys.stderr.write("examples are provided in the Git repo\n")
    sys.exit()


# Directory paths
DBDIR = config['Paths']['Databases']
TMPDIR = config['Paths']['Temp']
# paths for our databases
BACPATH  = os.path.join(DBDIR, "bac_giant_unique_species")
HOSTPATH = os.path.join(DBDIR, "human_masked")
CONPATH  = os.path.join(DBDIR, "contaminants")
PROTPATH = os.path.join(DBDIR, "proteins")
TAXPATH  = os.path.join(DBDIR, "taxonomy")
NUCLPATH = os.path.join(DBDIR, "nucleotides")
URVPATH = os.path.join(PROTPATH, "uniref_plus_virus") # uniref50 + viruses


# create directories
for DIR in [DBDIR, TMPDIR, BACPATH, HOSTPATH, CONPATH, PROTPATH, TAXPATH, NUCLPATH, URVPATH]:
    if not os.path.exists(DIR):
        os.mkdir(DIR)


# import the rules
include: "rules/99_download_dbs.smk"


## bowtie vs bbmap databases
inputs = []
if config['Options']['use_bowtie']:
    inputs = [
        expand(os.path.join(BACPATH, "bac_uniquespecies_giant.masked_Ns_removed.{n}.bt2l"), n=[1,2,3,4]),
        expand(os.path.join(HOSTPATH, "human_virus_masked.{n}.bt2l"), n=[1,2,3,4]),
        expand(os.path.join(CONPATH, "line_sine.{n}.bt2"), n=[1,2,3,4]),
        expand(os.path.join(CONPATH, "line_sine.rev.{m}.bt2"), m=[1,2])
    ]
else:
    inputs = [
        os.path.join(BACPATH, "ref"),
        os.path.join(HOSTPATH, "ref"),
        os.path.join(CONPATH, "line_sine.fasta")
    ]


rule all:
    input:
        # the database directories
        inputs,
        #os.path.join(PROTPATH, "uniprot_virus.faa"),
        os.path.join(TAXPATH, "uniprot_ncbi_mapping.dat"),
        os.path.join(NUCLPATH, "refseq_virus_nt_UniVec_masked/nt.fnaDB.dbtype"),
        os.path.join(NUCLPATH, "bac_virus_masked/nt.fnaDB.dbtype"),
        os.path.join(NUCLPATH, "refseq_virus_nt_UniVec_masked/nt.fnaDB.index"),
        os.path.join(NUCLPATH, "bac_virus_masked/nt.fnaDB.index"),
        os.path.join(TAXPATH, "taxonomizr_accessionTaxa.sql"),
        multiext(os.path.join(PROTPATH, "uniprot_virus_c99"), ".db_mapping", ".db_names.dmp", ".db_nodes.dmp", ".db_merged.dmp", ".db_delnodes.dmp"),
        multiext(os.path.join(URVPATH, "uniref50_virus"), ".db_mapping", ".db_names.dmp", ".db_nodes.dmp", ".db_merged.dmp", ".db_delnodes.dmp")

