#!/usr/bin/env python3
"""Assemble the supplementary sequencing tables (sample inventory, replicon
coverage, mutation catalogue) as a single Excel workbook plus flat TSVs."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent

CONDITION_ORDER = [
    "ancestral clone",
    "TriComm KAN no drug", "TriComm KAN 16xMIC AMP",
    "TriComm KAN+PVI no drug", "TriComm KAN+PVI 16xMIC AMP",
    "TriComm PVI no drug", "TriComm PVI 16xMIC AMP",
    "TriComm PVB no drug", "TriComm PVB 16xMIC AMP",
    "TriComm KAN 1xMIC COL", "TriComm KAN 2xMIC COL",
    "TriComm PN23 no drug", "TriComm PN23 1xMIC COL", "TriComm PN23 2xMIC COL",
    "TriComm KAN+PN23 2xMIC COL",
]

CLASS_LABEL = {
    "candidate": "de novo",
    "ancestral background": "present in ancestor",
    "masked repeat/prophage": "masked (cryptic prophage / repeat)",
    "absent replicon": "replicon absent from sample",
    "engineered knockout": "engineered auxotrophy knockout",
}


def sort_table(df: pd.DataFrame, *extra: str) -> pd.DataFrame:
    """Order rows by experiment (ancestral, PN23, pOXA48), vial, replicate and
    sample, then by any extra columns."""
    cols = ["Experiment", "Vial", "Replicate", "Sample", *extra]
    return df.sort_values(
        cols, key=lambda c: c.str.lower() if c.name == "Experiment" else c
    )


def change_string(row) -> str:
    kind = row["kind"]
    ref = row.get("ref_seq")
    alt = row.get("detail")
    if kind == "SNP":
        return f"{ref}{int(row['position'])}{alt}"
    if kind == "SUB":
        return f"{alt} bp substitution at {int(row['position'])}"
    if kind == "DEL":
        return f"{alt} bp deletion at {int(row['position'])}"
    if kind == "INS":
        return f"+{alt} at {int(row['position'])}"
    return f"{kind} at {int(row['position'])}"


def aa_change(row) -> str:
    if pd.isna(row.get("aa_ref_seq")) or pd.isna(row.get("aa_position")):
        return ""
    pos = str(row["aa_position"])
    try:
        pos = str(int(float(pos)))
    except ValueError:
        pass
    return f"{row['aa_ref_seq']}{pos}{row['aa_new_seq']}"


def main() -> None:
    ident = pd.read_csv(OUT / "sample_identity.tsv", sep="\t")
    mut = pd.read_csv(OUT / "mutations_classified.tsv", sep="\t", low_memory=False)

    # ------------------------------------------------------------- sample sheet
    samples = ident.rename(
        columns={
            "sample_id": "Sample",
            "experiment": "Experiment",
            "replicate": "Replicate",
            "vial": "Vial",
            "condition": "Condition",
            "declared_strain": "Expected strain",
            "sample_type": "Sample type",
            "timepoint_h": "Time (h)",
            "chromosome": "Chromosome coverage (x)",
            "genotype_strain": "Genotype from auxotrophy deletions",
            "deleted_aux_genes": "Deleted auxotrophy genes",
            "observed_marker": "Fluorescent marker from coverage",
            "msfGFP_copy_number": "pS2313M (msfGFP) / chromosome",
            "mScarlet_copy_number": "pS2313RS (mScarlet) / chromosome",
            "pOXA48_copy_number": "pOXA-48 / chromosome",
            "PN23_copy_number": "PN23 / chromosome",
            "genotype_check": "Auxotrophy genotype check",
            "marker_check": "Marker plasmid check",
        }
    )
    keep = [
        "Sample", "Experiment", "Replicate", "Vial", "Condition", "Sample type",
        "Expected strain", "Time (h)", "Chromosome coverage (x)",
        "Deleted auxotrophy genes", "Genotype from auxotrophy deletions",
        "Auxotrophy genotype check", "Fluorescent marker from coverage",
        "Marker plasmid check", "pS2313M (msfGFP) / chromosome",
        "pS2313RS (mScarlet) / chromosome", "pOXA-48 / chromosome",
        "PN23 / chromosome",
    ]
    samples = samples[keep].copy()
    for col in keep[-4:] + ["Chromosome coverage (x)"]:
        samples[col] = samples[col].astype(float).round(2)

    # The auxotrophy-deletion genotype is only reported for a population when
    # the deletion pattern matches a single strain (i.e. the population is
    # effectively clonal); a mixed community gives an ambiguous pattern that is
    # not a genotype call, so those cells are left empty.
    population = samples["Sample type"].eq("population")
    ambiguous = samples["Genotype from auxotrophy deletions"].eq("ambiguous")
    samples.loc[population & ambiguous,
                ["Deleted auxotrophy genes", "Genotype from auxotrophy deletions"]] = ""
    # Identity checks against the declared strain are only meaningful for the
    # evolved clones: leave them empty for ancestors and populations.
    unchecked = population | samples["Experiment"].eq("ancestral")
    samples.loc[unchecked, ["Auxotrophy genotype check",
                            "Fluorescent marker from coverage"]] = ""
    samples["Marker plasmid check"] = samples["Marker plasmid check"].replace(
        "not applicable", "N/A"
    )
    samples = sort_table(samples)

    # ----------------------------------------------------------- mutation sheet
    m = mut[mut["call_class"] != "absent replicon"].copy()
    m["Change"] = m.apply(change_string, axis=1)
    m["Amino acid change"] = m.apply(aa_change, axis=1)
    m["Class"] = m["call_class"].map(CLASS_LABEL)
    m["Detection"] = np.where(m["polymorphism_mode"], "polymorphism mode", "consensus")
    m["Frequency"] = pd.to_numeric(m["frequency"], errors="coerce")
    m.loc[m["Detection"].eq("consensus"), "Frequency"] = 1.0

    table = m.rename(
        columns={
            "sample_id": "Sample", "experiment": "Experiment", "replicate": "Replicate",
            "vial": "Vial", "condition": "Condition", "declared_strain": "Strain",
            "sample_type": "Sample type", "timepoint_h": "Time (h)",
            "replicon_name": "Replicon", "position": "Position", "kind": "Type",
            "gene_name": "Gene", "gene_position": "Gene position",
            "gene_product": "Product", "ancestors_carrying": "Ancestor carrying this call", "mutation_category": "Category",
        }
    )[
        ["Sample", "Experiment", "Replicate", "Vial", "Condition", "Sample type",
         "Strain", "Time (h)", "Replicon", "Position", "Type", "Change",
         "Gene", "Amino acid change", "Gene position", "Category", "Product",
         "Detection", "Frequency", "Class", "Ancestor carrying this call"]
    ]
    table = sort_table(table, "Replicon", "Position")

    # all de novo calls are kept, including sites seen only at low frequency in
    # many unrelated samples; the per-site sample counts are in table S4
    denovo = table[table["Class"] == "de novo"]

    # --------------------------------------------------- per-gene parallelism
    par = (
        denovo.groupby(["Gene", "Change", "Amino acid change", "Category"])
        .agg(samples=("Sample", "nunique"),
             consensus_samples=("Detection", lambda s: (s == "consensus").sum()),
             max_frequency=("Frequency", "max"),
             conditions=("Condition", lambda s: "; ".join(sorted(set(s.dropna())))))
        .reset_index()
        .sort_values(["samples", "max_frequency"], ascending=False)
    )

    with pd.ExcelWriter(OUT / "supplementary_sequencing_tables.xlsx") as xl:
        samples.to_excel(xl, sheet_name="S1 samples and coverage", index=False)
        table.to_excel(xl, sheet_name="S2 all mutation calls", index=False)
        denovo.to_excel(xl, sheet_name="S3 de novo mutations", index=False)
        par.to_excel(xl, sheet_name="S4 recurrent mutations", index=False)

    samples.to_csv(OUT / "table_S1_samples_coverage.tsv", sep="\t", index=False)
    table.to_csv(OUT / "table_S2_all_mutations.tsv", sep="\t", index=False)
    denovo.to_csv(OUT / "table_S3_denovo_mutations.tsv", sep="\t", index=False)
    par.to_csv(OUT / "table_S4_recurrent_mutations.tsv", sep="\t", index=False)

    print(f"samples: {len(samples)}")
    print(f"all calls (excl. absent replicons): {len(table)}")
    print(f"de novo candidate calls: {len(denovo)}")
    print(f"distinct de novo mutations: {len(par)}")
    print("\nde novo calls per class of sample:")
    print(denovo.groupby(["Experiment", "Sample type"]).size().to_string())


if __name__ == "__main__":
    main()
