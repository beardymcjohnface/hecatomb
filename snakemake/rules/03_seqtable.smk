"""

Snakefile specifically to run the R code.

Rob Edwards, Feb 2020
"""

rule merge_seq_table:
    """
    Reads the sequences and counts from each samples seqtable.txt file and converts to fasta format
    for the rest of the pipline. 
    """
    input:
        expand(os.path.join(QC, "counts", "{sample}_seqtable.txt"), sample=SAMPLES)
    output:
        fa = os.path.join(RESULTS,"seqtable.fasta")
    params:
        resultsdir = directory(RESULTS),
    benchmark:
        BENCHDIR + "/merge_seq_table.txt"
    run:
        out = open(output[0], 'w')
        seqId = 0
        for sample in SAMPLES:
            counts = open(os.path.join(QC, 'counts', f'{sample}_seqtable.txt'), 'r')
            # skip header
            line = counts.readline()
            for line in counts:
                id = str(seqId).zfill(8)
                seqId = seqId + 1
                l = line.split()
                out.write(f'>{id}\ts={sample}\tc={l[1]}\n{l[0]}\n')
            counts.close()
        out.close()
