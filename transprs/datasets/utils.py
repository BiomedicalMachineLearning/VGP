default_cnames = {
    # CHROMOSOME
    "CHR": "CHR",
    # CHR coordinate
    "BP": "BP",
    # RS NUMBER
    "SNP": "SNP",
    "MARKERNAME": "SNP",
    "SNPID": "SNP",
    "RS": "SNP",
    "RSID": "SNP",
    "RS_NUMBER": "SNP",
    "RS_NUMBERS": "SNP",
    # NUMBER OF STUDIES
    "NSTUDY": "NSTUDY",
    "N_STUDY": "NSTUDY",
    "NSTUDIES": "NSTUDY",
    "N_STUDIES": "NSTUDY",
    # P-VALUE
    "P": "P",
    "PVALUE": "P",
    "P_VALUE": "P",
    "PVAL": "P",
    "P_VAL": "P",
    "GC_PVALUE": "P",
    # ALLELE 1
    "A1": "A1",
    "ALLELE1": "A1",
    "ALLELE_1": "A1",
    "EFFECT_ALLELE": "A1",
    "REFERENCE_ALLELE": "A1",
    "INC_ALLELE": "A1",
    "EA": "A1",
    # ALLELE 2
    "A2": "A2",
    "ALLELE2": "A2",
    "ALLELE_2": "A2",
    "OTHER_ALLELE": "A2",
    "NON_EFFECT_ALLELE": "A2",
    "DEC_ALLELE": "A2",
    "NEA": "A2",
    # N
    "N": "N",
    "NCASE": "N_CAS",
    "CASES_N": "N_CAS",
    "N_CASE": "N_CAS",
    "N_CASES": "N_CAS",
    "N_CONTROLS": "N_CON",
    "N_CAS": "N_CAS",
    "N_CON": "N_CON",
    "N_CASE": "N_CAS",
    "NCONTROL": "N_CON",
    "CONTROLS_N": "N_CON",
    "N_CONTROL": "N_CON",
    "WEIGHT": "N",  # metal does this. possibly risky.
    # SIGNED STATISTICS
    "ZSCORE": "Z",
    "Z-SCORE": "Z",
    "GC_ZSCORE": "Z",
    "Z": "Z",
    "OR": "OR",
    "B": "BETA",
    "BETA": "BETA",
    "LOG_ODDS": "LOG_ODDS",
    "EFFECTS": "BETA",
    "EFFECT": "BETA",
    "SIGNED_SUMSTAT": "SIGNED_SUMSTAT",
    # INFO
    "INFO": "INFO",
    # MAF
    "EAF": "FRQ",
    "FRQ": "FRQ",
    "MAF": "FRQ",
    "FRQ_U": "FRQ",
    "F_U": "FRQ",
    # SE
    "SE": "SE",
}


def pre_reader(sumstat):
    clean_ignore = [clean_header(x) for x in sumstat.columns]

    cleaned_column = [default_cnames[x] for x in sumstat.columns]

    assert (
        "P" in cleaned_column
    ), "We cannot detect the P (p-value) column in the sumstat!"
    assert (
        "CHR" in cleaned_column
    ), "We cannot detect the CHR (Chromosome information) column in the sumstat!"
    assert (
        "SNP" in cleaned_column
    ), "We cannot detect the SNP (SNP information) column in the sumstat!"
    assert (
        "A1" in cleaned_column
    ), "We cannot detect the A1 (Allele 1) column in the sumstat!"
    assert (
        "A2" in cleaned_column
    ), "We cannot detect the A2 (Allele 2) column in the sumstat!"
    assert (
        "OR" or "BETA"
    ) in cleaned_column, (
        "We cannot detect the OR/BETA (Odd ratio/Beta) column in the sumstat!"
    )

    sumstat.columns = cleaned_column

    if "OR" in cleaned_column:
        default_order = [
            "CHR",
            "BP",
            "SNP",
            "A1",
            "A2",
            "N",
            "SE",
            "P",
            "OR",
            "INFO",
            "FRQ",
        ]
    else:
        default_order = [
            "CHR",
            "BP",
            "SNP",
            "A1",
            "A2",
            "N",
            "SE",
            "P",
            "BETA",
            "INFO",
            "FRQ",
        ]

    sumstat = sumstat.reindex(default_order, axis=1)

    return sumstat


def clean_header(header):
    """
    For cleaning file headers.
    - convert to uppercase
    - replace dashes '-' with underscores '_'
    - replace dots '.' (as in R) with underscores '_'
    - remove newlines ('\n')
    """
    return header.upper().replace("-", "_").replace(".", "_").replace("\n", "")


def complement(x):
    switch = {"A": "T", "T": "A", "G": "C", "C": "G"}
    if x in switch:
        return switch[x]
    else:
        return x
