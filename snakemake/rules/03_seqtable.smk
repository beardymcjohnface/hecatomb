"""

Snakefile specifically to run the R code.

Rob Edwards, Feb 2020
"""

rule merge_seq_table:
    """
    Merge seq counts
    """
    input:
        expand(os.path.join(QC, "counts", "{sample}_seqtable.txt"), sample=SAMPLES)
    output:
        os.path.join(RESULTS,"seqtable.fasta")
    params:
        resultsdir = directory(RESULTS),
    benchmark:
        "benchmarks/merge_seq_table.txt"
    run:
        out = open(output[0], 'w')
        seqId = 0
        for sample in SAMPLES:
            counts = open(os.path.join(QC, 'counts', f'{sample}_seqtable.txt'), 'r')
            for line in counts:
                id = str(seqId).zfill(8)
                seqId = seqId + 1
                l = line.split()
                out.write(f'>{id} s={sample} c={l[1]}\n{l[0]}\n')
            counts.close()
        out.close()
