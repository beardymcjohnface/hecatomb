rule map_shred_seq:
    """
    Map the virus shred seqs to the input host
    """
    input:
        ref = hostFasta,
        shred = SHRED
    output:
        temp(hostFasta + 'sam.gz')
    shell:
        """
        bbmap.sh ref={input.ref} in={input.shred} \
            outm={output} path=tmp/ \
            minid=0.90 maxindel=2 ow=t
        """


rule mask_host:
    """
    Mask the host
    """
    input:
        ref = hostFasta,
        sam = hostFasta + 'sam.gz'
    output:
        hostOutFasta
    shell:
        """
        bbmask.sh in={input.ref} out={output} entropy={entropy} sam={input.sam} ow=t
        """
