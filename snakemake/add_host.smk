
import os
import sys


# REQUIRED CONFIG
DBDIR = config['Databases']
hostFasta = config['HostFa']
hostName = config['HostName']
entropy = config['Entropy']

# OUTPUT HOST FASTA
hostOutFasta = os.path.join(DBDIR, 'hosts', hostName + '_virus_masked.fasta')


# IMPORT THE RULES
include: "rules/99_mask_host.smk"


# RUN
rule all:
    input:
        hostOutFasta
