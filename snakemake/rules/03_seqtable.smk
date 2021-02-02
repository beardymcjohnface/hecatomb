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
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        cts = os.path.join(RESULTS, 'seqCounts.tsv')
    params:
        resultsdir = directory(RESULTS),
    benchmark:
        "benchmarks/merge_seq_table.txt"
    run:
        outFa = open(output.fa, 'w')
        outCt = open(output.cts, 'w')
        seqId = 0
        for sample in SAMPLES:
            counts = open(os.path.join(QC, 'counts', f'{sample}_seqtable.txt'), 'r')
            for line in counts:
                id = str(seqId).zfill(8)
                seqId = seqId + 1
                l = line.split()
                outFa.write(f'>{id}\n{l[0]}\n')
                outCt.write(f'{id}\t{sample}\t{l[1]}\n')
            counts.close()
        outFa.close()
        outCts.close()
