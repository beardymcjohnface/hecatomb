
import os
import sys


# REQUIRED CONFIG
DBDIR = os.path.join(workflow.basedir, '../databases')
virShred = os.path.join(DBDIR, 'hosts', 'virus_shred.fasta.gz')
hostFasta = config['HostFa']
hostName = config['HostName']
entropy = config['Entropy']


# FOR SLURM
LOGDIR = 'logs'
if not os.path.exists(LOGDIR):
    os.mkdir(LOGDIR)


# OUTPUT HOST FASTA
hostOutFasta = os.path.join(DBDIR, 'hosts', f'{hostName}_virus_masked.fasta')


# IMPORT THE RULES
include: "rules/99_mask_host.smk"


# RUN
rule all:
    input:
        hostOutFasta
