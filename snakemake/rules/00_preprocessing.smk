"""
Snakemake rule file to preprocess Illumina sequence data in preperation for taxonomic assignment.

What is accomplished with this script?
    - Non-biological sequence removal (primers, adapters)

Additional Reading:
    - Hecatomb GitHub: https://github.com/shandley/hecatomb
    - Official Snakemake documentation: https://snakemake.readthedocs.io/en/stable/

Historical Notes:
- Updated Snakefile based on [contaminant_removal.sh](../base/contaminant_removal.sh)

Rob Edwards, Jan 2020
Updated: Scott Handley, Jan 2021
"""

# NOTE: bbtools uses "threads=auto" by default that typically uses all threads, so no need to specify. 
# -Xmx is used to specify the memory allocation for bbtools operations
# Set your -Xmx specifications in your configuration file 


rule remove_leftmost_primerB:
    """
    Step 01: Remove leftmost primer.
    Primer sequences used in the Handley lab are included (primerB.fa). If your lab uses other primers you will need to 
    place them in CONPATH (defined in the Snakefile) and change the file name from primerB.fa to your file name below.
    """
    input:
        r1 = os.path.join(READDIR, PATTERN_R1 + file_extension),
        r2 = os.path.join(READDIR, PATTERN_R2 + file_extension),
        primers = os.path.join(CONPATH, "primerB.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_01", f"{PATTERN_R1}.s1.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_01", f"{PATTERN_R2}.s1.out.fastq")),
        stats = os.path.join(STATS, "step_01", "{sample}.s1.stats.tsv")
    benchmark:
        os.path.join(BENCHDIR, "remove_leftmost_primerB.{sample}.txt")
    log:
        os.path.join(STDERR, "remove_leftmost_primerB.{sample}.log")
    resources:
        mem_mb=64000,
        cpus=16
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=11 ktrim=l restrictleft=20 \
            removeifeitherbad=f trimpolya=10 ordered=t rcomp=f ow=t \
            threads={resources.cpus} -Xmx{resources.mem_mb}m 2> {log}
        """


rule remove_3prime_contaminant:
    """
    Step 02: Remove 3' read through contaminant
    """
    input:
        r1 = os.path.join(TMPDIR, "step_01", f"{PATTERN_R1}.s1.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_01", f"{PATTERN_R2}.s1.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_02", f"{PATTERN_R1}.s2.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_02", f"{PATTERN_R2}.s2.out.fastq")),
        stats = os.path.join(STATS, "step_02", "{sample}.s2.stats.tsv")
    benchmark:
        os.path.join(BENCHDIR, "remove_3prime_contaminant.{sample}.txt")
    log:
        os.path.join(STDERR, "remove_3prime_contaminant.{sample}.log")
    resources:
        mem_mb=64000,
        cpus=16
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=11 ktrim=r removeifeitherbad=f ordered=t rcomp=f ow=t \
            threads={resources.cpus} -Xmx{resources.mem_mb}m 2> {log}
        """


rule remove_primer_free_adapter:
    """
    Step 03: Remove primer free adapter (both orientations)
    """
    input:
        r1 = os.path.join(TMPDIR, "step_02", f"{PATTERN_R1}.s2.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_02", f"{PATTERN_R2}.s2.out.fastq"),
        primers = os.path.join(CONPATH, "nebnext_adapters.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_03", f"{PATTERN_R1}.s3.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_03", f"{PATTERN_R2}.s3.out.fastq")),
        stats = os.path.join(STATS, "step_03", "{sample}.s3.stats.tsv")
    benchmark:
        os.path.join(BENCHDIR, "remove_primer_free_adapter.{sample}.txt")
    log:
        os.path.join(STDERR, "remove_primer_free_adapter.{sample}.log")
    resources:
        mem_mb=64000,
        cpus=16
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=10 ktrim=r removeifeitherbad=f ordered=t rcomp=t ow=t \
            threads={resources.cpus} -Xmx{resources.mem_mb}m 2> {log}
        """

rule remove_adapter_free_primer:
    """
    Step 04: Remove adapter free primer (both orientations)
    """
    input:
        r1 = os.path.join(TMPDIR, "step_03", f"{PATTERN_R1}.s3.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_03", f"{PATTERN_R2}.s3.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_04", f"{PATTERN_R1}.s4.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_04", f"{PATTERN_R2}.s4.out.fastq")),
        stats = os.path.join(STATS, "step_04", "{sample}.s4.stats.tsv")
    benchmark:
        os.path.join(BENCHDIR, "remove_adapter_free_primer.{sample}.txt")
    log:
        os.path.join(STDERR, "remove_adapter_free_primer.{sample}.log")
    resources:
        mem_mb=64000,
        cpus=16
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=0 removeifeitherbad=f ordered=t rcomp=t ow=t \
            threads={resources.cpus} -Xmx{resources.mem_mb}m 2> {log}
        """


rule remove_vector_contamination:
    """
    Step 05: Vector contamination removal (PhiX + NCBI UniVecDB)
    """
    input:
        r1 = os.path.join(TMPDIR, "step_04", f"{PATTERN_R1}.s4.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_04", f"{PATTERN_R2}.s4.out.fastq"),
        primers = os.path.join(CONPATH, "vector_contaminants.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_05", f"{PATTERN_R1}.s5.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_05", f"{PATTERN_R2}.s5.out.fastq")),
        stats = os.path.join(STATS, "step_05", "{sample}.s5.stats.tsv")
    benchmark:
        os.path.join(BENCHDIR, "remove_vector_contamination.{sample}.txt")
    log:
        os.path.join(STDERR, "remove_vector_contamination.{sample}.log")
    resources:
        mem_mb=64000,
        cpus=16
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=31 hammingdistance=1 ordered=t ow=t \
            threads={resources.cpus} -Xmx{resources.mem_mb}m 2> {log}
        """


rule remove_low_quality:
    """
    Step 06: Remove remaining low-quality bases and short reads
    """
    input:
        r1 = os.path.join(TMPDIR, "step_05", f"{PATTERN_R1}.s5.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_05", f"{PATTERN_R2}.s5.out.fastq")
    output:
        r1 = temp(os.path.join(TMPDIR, f"{PATTERN_R1}.clean.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, f"{PATTERN_R2}.clean.out.fastq")),
        stats = os.path.join(STATS, "step_06", "{sample}.s6.stats.tsv")
    benchmark:
        os.path.join(BENCHDIR, "remove_low_quality.{sample}.txt")
    log:
        os.path.join(STDERR, "remove_low_quality.{sample}.log")
    resources:
        mem_mb=64000,
        cpus=16
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            ordered=t qtrim=r maxns=2 \
            entropy={config[ENTROPY]} \
            trimq={config[QSCORE]} \
            minlength={config[MINLENGTH]} \
            threads={resources.cpus} -Xmx{resources.mem_mb}m 2> {log} 
        """


rule create_host_index:
    """
    Create the minimap2 index for mapping to the host; this will save a ton of time.
    """
    input:
        HOSTPATH,
        os.path.join(CONPATH, "line_sine.fasta")
    output:
        HOSTINDEX
    benchmark:
        os.path.join(BENCHDIR, "create_host_index.txt")
    log:
        os.path.join(STDERR, 'create_host_index.log')
    resources:
        mem_mb=64000,
        cpus=16
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -d {output} <(cat {input})"


rule host_removal_mapping:
    """
    Step 07a: Host removal. Must define host in config file (see Paths: Host: in config.yaml). Host should be masked of 
    viral sequence.
    Step 7 has several substeps (7a, 7b, 7c and 7d) to subset and fix read pairing issues
    If your reference is not available post an issue on GitHub requesting it to be added 
    (https://github.com/shandley/hecatomb)
    """
    input:
        r1 = os.path.join(TMPDIR, f"{PATTERN_R1}.clean.out.fastq"),
        r2 = os.path.join(TMPDIR, f"{PATTERN_R2}.clean.out.fastq"),
        host = HOSTINDEX
    output:
        os.path.join(QC, "HOST_REMOVED", f"{PATTERN_R1}.all.fastq")
    benchmark:
        os.path.join(BENCHDIR, "host_removal_mapping.{sample}.txt")
    log:
        mm=os.path.join(STDERR, "host_removal_mapping.{sample}.minimap.log"),
        sv=os.path.join(STDERR, "host_removal_mapping.{sample}.samtoolsView.log"),
        fq=os.path.join(STDERR, "host_removal_mapping.{sample}.samtoolsFastq.log")
    resources:
        mem_mb=64000,
        cpus=16
    conda:
        "../envs/minimap2.yaml"
    shell: # todo, fix this bandaid solution for paired fastq files
        """
        minimap2 -ax sr -t {resources.cpus} --secondary=no \
            {input.host} {input.r1} {input.r2} 2> {log.mm} \
            | samtools view -f 4 -h 2> {log.sv} \
            | samtools fastq -NO > {output} 2> {log.fq}
        """


# rule nonhost_read_repair:
#     """
#     Step 07c: Parse R1/R2 singletons (if singletons at all)
#     """
#     input:
#         singletons = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".unmapped.singletons.fastq")
#     output:
#         r1 = temp(os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".singletons.fastq")),
#         r2 = temp(os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".singletons.fastq"))
#     benchmark:
#         BENCHDIR + "/preprocessing/step_07/nonhost_read_repair_{sample}.txt"
#     log:
#         STDERR + "/host_removal/{sample}.nonhost_read_repair.log"
#     resources:
#         mem_mb=8000,
#         cpus=1
#     conda:
#         "../envs/bbmap.yaml"
#     shell:
#         """
#         reformat.sh in={input.singletons} out={output.r1} out2={output.r2} \
#             -Xmx{resources.mem_mb}m 2> {log}
#         """
#
#
# rule nonhost_read_combine:
#     """
#     Step 07d: Combine R1+R1_singletons and R2+R2_singletons
#     """
#     input:
#         r1 = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".unmapped.fastq"),
#         r2 = os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".unmapped.fastq"),
#         r1s = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".singletons.fastq"),
#         r2s = os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".singletons.fastq")
#     output:
#         r1 = temp(os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".all.fastq")),
#         r2 = temp(os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".all.fastq"))
#     benchmark:
#         BENCHDIR + "/preprocessing/step_07/nonhost_read_combine_{sample}.txt"
#     # resources:
#     #     mem_mb=100000,
#     #     cpus=64
#     shell:
#         """
#         cat {input.r1} {input.r1s} > {output.r1};
#         cat {input.r2} {input.r2s} > {output.r2}
#         """


rule remove_exact_dups:
    """
    Step 08: Remove exact duplicates
    """
    input:
        os.path.join(QC, "HOST_REMOVED", f"{PATTERN_R1}.all.fastq")
    output:
        temp(os.path.join(QC, "CLUSTERED", f"{PATTERN_R1}.deduped.out.fastq"))
    benchmark:
        os.path.join(BENCHDIR, "remove_exact_dups.{sample}.txt")
    log:
        os.path.join(STDERR, "remove_exact_dups.{sample}.log")
    resources:
        mem_mb=64000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        dedupe.sh in={input} out={output} \
            ac=f ow=t \
            threads={resources.cpus} -Xmx{resources.mem_mb}m 2> {log}
        """


rule cluster_similar_sequences:
    """
    Step 09: Cluster similar sequences at CLUSTERID in config.yaml. Default: 97% identity
    """
    input:
        os.path.join(QC, "CLUSTERED", f"{PATTERN_R1}.deduped.out.fastq")
    output:
        temp(os.path.join(QC, "CLUSTERED", "LINCLUST", f"{PATTERN_R1}_rep_seq.fasta")),
        temp(os.path.join(QC, "CLUSTERED", "LINCLUST", f"{PATTERN_R1}_cluster.tsv")),
        temp(os.path.join(QC, "CLUSTERED", "LINCLUST", f"{PATTERN_R1}_all_seqs.fasta"))
    params:
        respath=os.path.join(QC, "CLUSTERED", "LINCLUST"),
        tmppath=os.path.join(QC, "CLUSTERED", "LINCLUST", "{sample}_TMP"),
        prefix=PATTERN_R1
    benchmark:
        os.path.join(BENCHDIR, "cluster_similar_sequences.{sample}.txt")
    log:
        os.path.join(STDERR, "cluster_similar_sequences.{sample}.log")
    resources:
        mem_mb=64000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """ 
        mmseqs easy-linclust {input} {params.respath}/{params.prefix} {params.tmppath} \
            --kmer-per-seq-scale 0.3 \
            -c 0.95 --cov-mode 1 --threads {resources.cpus} &>> {log}
        """


rule create_individual_seqtables:
    """
    Step 10: Create individual seqtables. A seqtable is a count table with each feature (sequence) as a row, each column 
    as a sample and each cell the counts of each sequence per sample
    """
    input:
        seqs=os.path.join(QC, "CLUSTERED", "LINCLUST", f"{PATTERN_R1}_rep_seq.fasta"),
        counts=os.path.join(QC, "CLUSTERED", "LINCLUST", f"{PATTERN_R1}_cluster.tsv")
    output:
        seqs=temp(os.path.join(QC, "CLUSTERED", "LINCLUST", f"{PATTERN_R1}.seqs")),
        counts=temp(os.path.join(QC, "CLUSTERED", "LINCLUST", f"{PATTERN_R1}.counts")),
        seqtable=temp(os.path.join(QC, "CLUSTERED", "LINCLUST", f"{PATTERN_R1}.seqtable"))
    benchmark:
        os.path.join(BENCHDIR, "create_individual_seqtables.{sample}.txt")
    resources:
        mem_mb=64000,
        cpus=16
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit sort {input.seqs} --quiet -j {resources.cpus} -w 5000 -t dna \
            | seqkit fx2tab -w 5000 -t dna \
            | sed 's/\\t\\+$//' \
            | cut -f2,3 \
            | sed '1i sequence' > {output.seqs};
        cut -f1 {input.counts} \
            | sort \
            | uniq -c \
            | awk -F ' ' '{{print$2"\\t"$1}}' \
            | cut -f2 \
            | sed "1i {wildcards.sample}" > {output.counts};
        paste {output.seqs} {output.counts} > {output.seqtable};
        """


# rule merge_individual_seqtables:
#     """
#     Step 11: Merge individual sequence tables into combined seqtable
#     """
#     input:
#         files = expand(os.path.join(QC, "CLUSTERED", "LINCLUST", PATTERN_R1 + ".seqtable"), sample=SAMPLES)
#     output:
#         seqtable = os.path.join(RESULTS, "seqtable_all.tsv"),
#         tab2fx = temporary(os.path.join(RESULTS, "seqtable.tab2fx"))
#     params:
#         resultsdir = directory(RESULTS),
#     benchmark:
#         BENCHDIR + "/preprocessing/step_10/merge_seq_table.txt"
#     # resources:
#     #     mem_mb=100000,
#     #     cpus=64
#     params:
#         resultsdir = directory(RESULTS),
#     conda:
#         "../envs/R.yaml"
#     script:
#         "../scripts/seqtable_merge.R"


# rule convert_seqtable_tab_2_fasta:
#     """
#     Step 12: Convert tabular seqtable output to fasta
#     """
#     input:
#         os.path.join(RESULTS, "seqtable.tab2fx")
#     output:
#         os.path.join(RESULTS, "seqtable.fasta")
#     benchmark:
#         BENCHDIR + "/preprocessing/step_10/convert_seqtable_tab_2_fasta.txt"
#     resources:
#         mem_mb=64000,
#         cpus=16
#     conda:
#         "../envs/seqkit.yaml"
#     shell:
#         """
#         seqkit tab2fx {input} -j {resources.cpus} -w 5000 -t dna -o {output}
#         """


rule merge_seq_table:
    """
    Reads the sequences and counts from each samples seqtable text file and converts to fasta format
    for the rest of the pipline. 
    """
    input:
        expand(os.path.join(QC, "CLUSTERED", "LINCLUST", "{sample}_R1.seqtable"), sample=SAMPLES)
    output:
        fa = os.path.join(RESULTS,"seqtable.fasta")
    params:
        resultsdir = directory(RESULTS),
    benchmark:
        os.path.join(BENCHDIR, "merge_seq_table.txt")
    run:
        out = open(output[0], 'w')
        seqId = 0
        for sample in SAMPLES:
            counts = open(os.path.join(QC, "CLUSTERED", "LINCLUST", f"{sample}_R1.seqtable"), 'r')
            line = counts.readline() # skip header
            for line in counts:
                id = 's' + str(seqId).zfill(8)
                seqId = seqId + 1
                l = line.split()
                out.write(f'>{id}\ts={sample}\tc={l[1]}\n{l[0]}\n')
            counts.close()
        out.close()


rule create_seqtable_index:
    """
    Step 13: Index seqtable.fasta for rapid samtools access in later steps
    """
    input:
        os.path.join(RESULTS, "seqtable.fasta")
    output:
        os.path.join(RESULTS, "seqtable.faidx")
    conda:
        "../envs/samtools.yaml"
    log:
        os.path.join(STDERR, 'create_seqtable_index.log')
    resources:
        mem_mb=8000
    shell:
        """
        samtools faidx {input} -o {output}
        """


rule calculate_seqtable_sequence_properties:
    """
    Step 14: Calculate additional sequence properties (ie. GC-content) per sequence
    """
    input:
        os.path.join(RESULTS, "seqtable.fasta")
    output:
        os.path.join(RESULTS, "seqtable_properties.tsv")
    benchmark:
        os.path.join(BENCHDIR, "calculate_seqtable_sequence_properties.txt")
    log:
        os.path.join(STDERR, 'calculate_seqtable_sequence_properties.log')
    resources:
        mem_mb=64000,
        cpus=16
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit fx2tab -j {resources.cpus} --gc -H {input} | cut -f1,4 > {output}
        """

