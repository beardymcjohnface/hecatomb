"""
concatenate_results.sh

take all the results and make some cool data sets from them!
"""


rule concatenate_results_first:
    input:
        os.path.join(RESULTS, "viruses_tax_table.tsv"),
        os.path.join(RESULTS, "phage_tax_table.tsv"),
        os.path.join(RESULTS, "aa.aln.m8"),
        os.path.join(RESULTS, "nt.aln.m8"),
        "family_reads"

rule jive_aa_annotation:
    input:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table.tsv")
    output:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table_edited.tsv")
    shell:
        """
        sed 's/uc_//g' {input} > {output}
        """

rule jive_nt_annotation:
    input:
        os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_checked_lineage.tsv")
    output:
        os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_checked_lineage_edited.tsv")
    shell:
        """
        tail -n +2 {input} | sed 's/NA/unknown/g; s/uc_//g' > {output}
        """

# this is the name of the rule in concatenate_results.sh and it is too cute to change
rule happily_marry_outputs:
    input:
        aa = os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table_edited.tsv"),
        nt = os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_checked_lineage_edited.tsv")
    output:
        temporary(os.path.join(RESULTS, "viruses_tax_table_tmp.tsv"))
    shell:
        """
        cat {input.aa} {input.nt} | sort -n -k 1 > {output}
        """


rule add_crown_to_marriage:
    # not really, its a title
    input:
        tbl = os.path.join(RESULTS, "viruses_tax_table_tmp.tsv"),
        tax = os.path.join(TAXPATH, '2020_07_27_Viral_Baltimore_full_classification_table_ICTV2019.txt')
    output:
        os.path.join(RESULTS, "viruses_tax_table.tsv")
    run:
        # read in taxon classifications
        taxon = open(input.tax,'r')
        line = taxon.readline()
        taxonClass = {}
        balt = {}
        for line in taxon:
            l = line.rstrip().split('\t')
            # attach baltimore class and group to Genus
            try:
                balt[l[4]]
            except KeyError:
                balt[l[4]] = [l[7],l[8]]
            # attach child taxon to parent designations e.g. key=species, value=genus
            for i in range(5, 0, -1):
                if l[i] == '' or l[i-1] == '':
                    continue
                try:
                    taxonClass[l[i]]
                except KeyError:
                    taxonClass[l[i]] = l[i-1]
        taxon.close()
        # parse the taxon results and fill in any blanks
        table = open(input.tbl, 'r')
        out = open(output[0],'w')
        out.write('id\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tBaltimore\tBaltimoreGroup')
        for line in table:
            l = line.rstrip().split('\t')
            for i in range(7,2,-1):
                if l[i-1] == 'unknown' and l[i] != 'unknown':
                    try:
                        taxonClass[l[i]]
                        l[i-1] = taxonClass[l[i]]
                    except KeyError:
                        pass
            try:
                l.append(balt[l[6]])
            except KeyError:
                l.append(['',''])
            out.write('\t'.join(l))
            out.write('\n')
        table.close()
        out.close()


rule fix_phage_names:
    input:
        tsv = os.path.join(AA_OUT, "phage_table.tsv")
    output:
        temporary(os.path.join(RESULTS, "phage_tax_table_temp.tsv"))
    shell:
        """
        sed 's/uc_//g' {input} > {output}
        """

rule fix_phage_cols:
    input:
        os.path.join(RESULTS, "phage_tax_table_temp.tsv")
    output:
        temporary(os.path.join(RESULTS, "phage_tax_table_temp2.tsv"))
    shell:
        """
        cut -f1-8 {input} > {output}
        """
            
rule add_phage_title:
    input:
        os.path.join(RESULTS, "phage_tax_table_temp2.tsv")
    output:
        os.path.join(RESULTS, "phage_tax_table.tsv")
    shell:
        """
        sed -e '1iid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' {input} > {output}
        """

rule add_aa_tax_header:
    input:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.m8")
    output:
        os.path.join(RESULTS, "aa.aln.m8")
    shell:
        """
        sed -e '1iquery\ttarget\tpercent_id\talignment_length\tnum_mismatches\tnumber_gaps\tstart_query\tend_query\tstart_target\tend_target\te_value\tbit_score' \
                {input} > {output}
        """


rule add_nt_tax_header:
    input:
        os.path.join(NT_OUT, "resultDB.firsthit.m8")
    output:
        os.path.join(RESULTS, "nt.aln.m8")
    shell:
        """
        sed -e '1iquery\ttarget\tpercent_id\talignment_length\tnum_mismatches\tnumber_gaps\tstart_query\tend_query\tstart_target\tend_target\te_value\tbit_score' \
                {input} > {output}
        """

"""
I am not very happy with this solution and looking for some alternative ideas
Basically this is very linear, and I can't figure out how to make it non-linear!
"""

rule make_fasta:
    input:
        os.path.join(RESULTS, "viruses_tax_table.tsv")
    output:
        directory("family_reads")
    shell:
        """
        mkdir -p {output} && 
        for FAM in $(tail -n +2 {input} | cut -f6 | awk '!s[$0]++');
        do 
            for TID in $(grep $FAM {input} | cut -f 1);
            do
                grep -A1 -Fw $TID results/seqtable.fasta >> {output}/$FAM.fasta
            done
        done
        """




