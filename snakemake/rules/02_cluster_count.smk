"""

Snakefile based on [cluster_count.sh](../base/cluster_count.sh)

Rob Edwards, Jan 2020
"""

import os
import sys

rule remove_exact_dups:
    """
    Step 1: Remove exact duplicates
    """
    input:
        os.path.join(QC, "step_9", PATTERN_R1 + ".viral_amb.fastq")
    output:
        temp(os.path.join(QC, "step_10", PATTERN_R1 + ".s9.deduped.out.fastq"))
    benchmark:
        BENCHDIR + "/remove_exact_dups_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        dedupe.sh in={input} \
                out={output} \
                ac=f  ow=t \
                -Xmx{resources.mem_mb}m
        """

rule deduplicate:
    """
    Step 2: Dereplicate
    """
    input:
        os.path.join(QC, "step_10", PATTERN_R1 + ".s9.deduped.out.fastq")
    output:
        fa = temp(os.path.join(QC, "step_11", PATTERN_R1 + ".best.fasta")),
        stats = temp(os.path.join(QC, "step_11", "{sample}_stats.txt"))
    benchmark:
        BENCHDIR + "/deduplicate_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        dedupe.sh in={input} \
            csf={output.stats} out={output.fa} \
            ow=t s=4 rnc=t pbr=t \
            -Xmx{resources.mem_mb}m
        """

rule extract_seq_counts:
    """
    Step 3: Make each sequences one line long
    """
    input:
        os.path.join(QC, "step_11", PATTERN_R1 + ".best.fasta")
    output:
        temp(os.path.join(QC, "step_12", PATTERN_R1 + ".reformated.fasta"))
    benchmark:
        BENCHDIR + "/extract_seq_counts_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        reformat.sh in={input} out={output} \
            deleteinput=t fastawrap=0 \
            ow=t \
            -Xmx{resources.mem_mb}m
        """

rule extract_counts:
    """
    put sequences into a text file
    """
    input:
        os.path.join(QC, "step_12", PATTERN_R1 + ".reformated.fasta")
    output:
        temp(os.path.join(QC, "counts", "{sample}_seqs.txt"))
    shell:
        """
        grep -v '>' {input} | sed '1i sequence' > {output}
        """

# rule extract_counts_ids:
#     """
#     put the sequence IDs into a list text file. This file is not used for anything.
#     """
#     input:
#         os.path.join(QC, "step_12", PATTERN_R1 + ".reformated.fasta")
#     output:
#         os.path.join(QC, "counts", "{sample}_contig_ids.txt")
#     shell:
#         """
#         grep '>' {input} | sed 's|>Cluster_||' | awk -F "," '{{ print$1 }}' | sort -n | sed '1i contig_ids' > {output}
#         """

rule exract_count_stats:
    """
    put the sequence counts into a list text file
    """
    input:
        os.path.join(QC, "step_11", "{sample}_stats.txt")
    output:
        temp(os.path.join(QC, "counts", "{sample}_counts.txt"))
    params:
        s = "{sample}"
    shell:
        """
        cut -f 2 {input} | sed "1s/size/{params.s}/" > {output}
        """

rule create_seq_table:
    """
    combine the sequences list and the counts list.
    """
    input:
        seq = os.path.join(QC, "counts", "{sample}_seqs.txt"),
        cnt = os.path.join(QC, "counts", "{sample}_counts.txt")
    output:
        os.path.join(QC, "counts", "{sample}_seqtable.txt")
    shell:
        """
        paste {input.seq} {input.cnt} > {output}
        """

