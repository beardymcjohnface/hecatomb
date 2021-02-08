

rule long_sample_table:
    """
    Make a long table with everything
    """
    input:
        tx = os.path.join(RESULTS, "viruses_tax_table.tsv"),
        fa = os.path.join(RESULTS,"seqtable.fasta"),
        nt = os.path.join(NT_OUT, "resultDB.firsthit.m8"),
        aa = os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.m8")
    output:
        os.path.join(RESULTS, 'big_table.tsv')
    resources:
        mem_mb=16000
    run:
        def slurpAli(fh):
            d = {}
            for line in fh:
                l = line.split()
                try:
                    d[l[0]]
                except KeyError:
                    d[l[0]] = [l[2], l[3], l[10]]
            return d
        # read nt alignments
        nt_fh = open(input.nt, 'r')
        nt_ali = slurpAli(nt_fh)
        nt_fh.close()
        # read aa alignments
        aa_fh = open(input.aa,'r')
        aa_ali = slurpAli(aa_fh)
        aa_fh.close()
        # read in samples and counts for each seq id
        counts = {}
        fa_fh = open(input.fa, 'r')
        for line in fa_fh:
            if line.startswith('>'):
                l = line.split()
                l[0] = l[0].replace('>','')
                l[1] = l[1].replace('s=','')
                l[2] = l[2].replace('c=','')
                try:
                    counts[l[0]]
                except KeyError:
                    counts[l[0]] = [l[1],l[2]]
        fa_fh.close()
        # parse and expand the virus tax table
        tx_fh = open(input.tx, 'r')
        line = tx_fh.readline()
        out = open(output[0],'w')
        out.write('id\tsample\tcount\tnt_aa\tali_perc\tali_len\tali_eval\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tBaltimore\tBaltimoreGroup\n')
        for line in tx_fh:
            l = line.rstrip().split('\t')
            try:
                l[1:1] = counts[l[0]]
            except KeyError:
                sys.stderr.write(f'seq ID {l[0]} not in {input.fa}???')
                exit(1)
            try:
                l.insert(3,aa_ali[l[0]] * 3)
                l.insert(3,'aa')
            except KeyError:
                try:
                    l.insert(3,nt_ali[l[0]])
                    l.insert(3,'nt')
                except KeyError:
                    sys.stderr.write(f'seq ID {l[0]} not in either {input.aa} or {input.nt}???')
                    exit(1)
            out.write('\t'.join(l))
            out.write('\n')
        out.close()
