"""
Scott Handley, Jan 2021
Updated by Michael Roach, Feb 2021
"""


rule assembly_kmer_normalization:
    """
    Step 15: Kmer normalization. Data reduction for assembly improvement
    """
    input:
        r1=os.path.join(QC,"HOST_REMOVED",PATTERN_R1 + ".unmapped.fastq"),
        r2=os.path.join(QC,"HOST_REMOVED",PATTERN_R2 + ".unmapped.fastq"),
        r1s=os.path.join(QC,"HOST_REMOVED",PATTERN_R1 + ".singletons.fastq"),
        r2s=os.path.join(QC,"HOST_REMOVED",PATTERN_R2 + ".singletons.fastq")
    output:
        r1_norm=temp(os.path.join(ASSEMBLY,PATTERN_R1 + ".norm.fastq")),
        r2_norm=temp(os.path.join(ASSEMBLY,PATTERN_R2 + ".norm.fastq"))
    benchmark:
        BENCHDIR + "/preprocessing/assembly/kmer_normalization_{sample}.txt"
    log:
        STDERR + "/assembly/{sample}.kmer_normalization.log"
    resources:
        mem_mb=64000,
        cpus=16
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbnorm.sh in={input.r1} in2={input.r2} \
            extra={input.r1s},{input.r2s} \
            out={output.r1_norm} out2={output.r2_norm} \
            target=100 ow=t \
            threads={resources.cpus} -Xmx{resources.mem_mb}m 2> {log}
    """


rule individual_sample_assembly:
    """
    Step 16: Individual sample assemblies
    """
    input:
        r1_norm=os.path.join(ASSEMBLY,PATTERN_R1 + ".norm.fastq"),
        r2_norm=os.path.join(ASSEMBLY,PATTERN_R2 + ".norm.fastq"),
        r1s=os.path.join(QC,"HOST_REMOVED",PATTERN_R1 + ".singletons.fastq"),
        r2s=os.path.join(QC,"HOST_REMOVED",PATTERN_R2 + ".singletons.fastq")
    output:
        contigs=os.path.join(ASSEMBLY,'{sample}',"{sample}.contigs.fa")
    params:
        mh_dir=directory(os.path.join(ASSEMBLY,'{sample}')),
        contig_dic=directory(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY"))
    benchmark:
        BENCHDIR + "/preprocessing/assembly/megahit_{sample}.txt"
    log:
        STDERR + "/assembly/{sample}.megahit.log"
    resources:
        mem_mb=64000,
        cpus=16
    conda:
        "../envs/megahit.yaml"
    shell:
        """
        rmdir {params.mh_dir};
        megahit -1 {input.r1_norm} -2 {input.r2_norm} -r {input.r1s},{input.r2s} \
            -o {params.mh_dir} --out-prefix {wildcards.sample} \
            --k-min 45 --k-max 225 --k-step 26 --min-count \
            -t {resources.cpus} &>> {log}      
        """


rule concatenate_contigs:
    """
    Step 17: Concatenate individual assembly outputs (contigs) into a single file
    """
    input:
        lambda wildcards: expand(os.path.join(ASSEMBLY,'{sample}','{sample}' + ".contigs.fa"),sample=SAMPLES)
    output:
        os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","all_megahit_contigs.fasta")
    shell:
        "cat {input} > {output}"


rule contig_reformating_and_stats:
    """
    Step 18: Remove short contigs (Default: 500). Defined in config[CONTIG_SIZE_THRESH]
    """
    input:
        os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","all_megahit_contigs.fasta")
    output:
        rename=temp(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","all_megahit_contigs.renamed.fasta")),
        size=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","all_megahit_contigs_size_selected.fasta"),
        stats=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","all_megahit_contigs.stats"),
        sketch=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","all_megahit_contigs_size_selected.sketch")
    benchmark:
        BENCHDIR + "/preprocessing/contig_reformating.txt"
    log:
        STDERR + "/assembly/contig_reformating.log"
    resources:
        mem_mb=8000
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        rename.sh in={input} out={output.rename} \
            ow=t \
            -Xmx{resources.mem_mb}m;
        reformat.sh in={output.rename} out={output.size} \
            ml={config[MINLENGTH]} \
            ow=t \
            -Xmx{resources.mem_mb}m;
        statswrapper.sh in={input} out={output.stats} \
            format=2 \
            ow=t;
        sendsketch.sh in={output.size} out={output.sketch} \
            address=nt mode=sequence format=3 \
            ow=t \
            -Xmx{resources.mem_mb}m 2> {log};
        """


rule population_assembly:
    """
    Step 19: Create 'contig dictionary' of all unique contigs present in the study (population assembly)
    """
    input:
        os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","all_megahit_contigs_size_selected.fasta")
    output:
        assembly=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","assembly.fasta"),
        stats=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","contig_dictionary.stats"),
        sketch=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","contig_dictionary.sketch")
    params:
        flye_out=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE")
    benchmark:
        BENCHDIR + "/assembly/population_assembly.txt"
    log:
        STDERR + "/assembly/population_assembly.log"
    resources:
        mem_mb=64000,
        cpus=16
    conda:
        "../envs/metaflye.yaml"
    shell:
        """
        flye --subassemblies {input} -t {resources.cpus} --plasmids -o {params.flye_out} -g 1g;
        statswrapper.sh in={output.assembly} out={output.stats} \
            format=2 ow=t;
        sendsketch.sh in={output.assembly} out={output.sketch} \
            address=nt mode=sequence format=3 ow=t \
            -Xmx{resources.mem_mb}m 2> {log};
        """


rule coverage_calculations:
    """
    Step 20a: Calculate contig coverage and extract unmapped reads
    """
    input:
        r1=os.path.join(QC,"HOST_REMOVED",PATTERN_R1 + ".all.fastq"),
        r2=os.path.join(QC,"HOST_REMOVED",PATTERN_R2 + ".all.fastq"),
        ref=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","assembly.fasta")
    output:
        sam=temp(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",'{sample}' + ".aln.sam.gz")),
        unmap=temp(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",'{sample}' + ".umapped.fastq")),
        covstats=temp(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",'{sample}' + ".cov_stats")),
        rpkm=temp(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",'{sample}' + ".rpkm")),
        statsfile=temp(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",'{sample}' + ".statsfile")),
        scafstats=temp(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",'{sample}' + ".scafstats"))
    benchmark:
        BENCHDIR + "/assembly/coverage_calculations_{sample}.txt"
    log:
        STDERR + "/assembly/coverage_calculations_{sample}.log"
    resources:
        mem_mb=64000,
        cpus=16
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbmap.sh ref={input.ref} in={input.r1} in2={input.r2} \
            nodisk \
            out={output.sam} \
            outu={output.unmap} \
            ambiguous=random \
            physcov=t \
            covstats={output.covstats} \
            rpkm={output.rpkm} \
            statsfile={output.statsfile} \
            scafstats={output.scafstats} \
            maxindel=100 minid=90 \
            ow=t 2> \
            threads={resources.cpus} -Xmx{resources.mem_mb}m 2> {log};
        """


rule create_contig_count_table:
    """
    Step 20b: Filter low coverage contigs
    """
    input:
        rpkm=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",'{sample}' + ".rpkm"),
        covstats=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",'{sample}' + ".cov_stats")
    output:
        counts_tmp=temp(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",'{sample}' + "_counts.tmp")),
        TPM_tmp=temp(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",'{sample}' + "_TPM.tmp")),
        TPM=temp(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",'{sample}' + "_TPM")),
        TPM_final=temp(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",'{sample}' + "_TPM.final")),
        cov_temp=temp(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",'{sample}' + "_cov.tmp")),
        count_tbl=temp(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",'{sample}' + "_contig_counts.tsv"))
    shell:
        """
        ## TPM Calculator
        # Prepare table & calculate RPK
        tail -n+6 {input.rpkm} > {output.counts_tmp} | \
        cut -f1,2,5,6,8 | \
        awk 'BEGIN{{ FS=OFS="\t" }} {{ print $0, $3/($2/1000) }}' > {output.counts_tmp}

        # Calculate size factor
        sizef=$(awk 'BEGIN{{ total=0 }} {{ total=total+$6 }} END{{ printf total }}' {output.counts_tmp});

        # Calculate TPM
        awk -v awkvar="$sizef" 'BEGIN{{ FS=OFS="\t" }} {{ print $0, $6/awkvar }}' < {output.counts_tmp} > {output.TPM_tmp};

        # Add sample name
        awk -v awkvar="{wildcards.sample}" 'BEGIN{{FS=OFS="\t"}} {{ print awkvar, $0 }}' < {output.TPM_tmp} > {output.TPM};

        # Remove RPK
        cut -f1-6,8 {output.TPM} > {output.TPM_final};

        ## Coverage stats modifications
        tail -n+2 {input.covstats} | cut -f2,4,5,6,9 > {output.cov_temp};

        ## Combine tables
        paste {output.TPM_final} {output.cov_temp} > {output.count_tbl};

        """


rule concatentate_contig_count_tables:
    """
    Rule 20c: Concatenate contig count tables
    Note: this is done as a separate rule due to how snakemake handles i/o files. It does not work well in Rule 20b as the i/o PATTERNS are different.
    """
    input:
        lambda wildcards: expand(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING",'{sample}' + "_contig_counts.tsv"),
            sample=SAMPLES)
    output:
        os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING","contig_count_table.tsv")
    shell:
        """
        cat {input} > {output};
        sed -i '1i sample_id\tcontig_id\tlength\treads\tRPKM\tFPKM\tTPM\tavg_fold_cov\tcontig_GC\tcov_perc\tcov_bases\tmedian_fold_cov' {output};
        """
