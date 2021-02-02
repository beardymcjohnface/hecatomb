

rule long_sample_table:
    """
    Make a long table per sample
    """
    input:
        balt = os.path.join(TAXPATH, '2020_07_27_Viral_Baltimore_full_classification_table_ICTV2019.txt'),
        idt = os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table_edited.tsv"),
        nt = os.path.join(NT_OUT, "resultDB.firsthit.m8"),
        aa = os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.m8")
    output: