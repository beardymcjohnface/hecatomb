rule map_shred_seq:
    """
    Map the virus shred seqs to the input host
    """
    input:
        ref = hostFasta,
        shred = virShred
    output:
        temp(f'{hostFasta}.sam.gz')
    conda:
        "../envs/bbmap.yaml"
    resources:
        mem_mb=64000,
        cpus=8
    shell:
        """
        bbmap.sh ref={input.ref} in={input.shred} \
            outm={output} path=tmp/ \
            minid=0.90 maxindel=2 ow=t \
            threads={resources.cpus} -Xmx{resources.mem_mb}m
        """


rule mask_host:
    """
    Mask the host
    """
    input:
        ref = hostFasta,
        sam = f'{hostFasta}.sam.gz'
    output:
        hostOutFasta
    conda:
        "../envs/bbmap.yaml"
    resources:
        mem_mb=64000,
        cpus=8
    shell:
        """
        bbmask.sh in={input.ref} out={output} \
            entropy={entropy} sam={input.sam} ow=t \
            threads={resources.cpus} -Xmx{resources.mem_mb}m
        """
