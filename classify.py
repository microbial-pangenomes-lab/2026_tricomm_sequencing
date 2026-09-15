#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parent

AUX_GENES = {
    "lysA": (2_970_996, 2_972_258),
    "pheA": (2_731_104, 2_732_264),
    "metA": (4_204_208, 4_205_137),
}
EXPECTED_KO = {"PL": {"pheA", "lysA"}, "PM": {"pheA", "metA"}, "LM": {"lysA", "metA"}}

# Cryptic-prophage / repeat loci in BW25113 that produce dense clusters of
# low-confidence calls. Masked from the adaptive-mutation analysis.
MASKED_REGIONS = [
    ("DLP12 prophage (ylcG-ompT)", 563_000, 586_000),
    ("e14 prophage (ymfH-ymfR)", 1_195_000, 1_202_500),
    ("Rac prophage (stfR)", 1_422_500, 1_424_500),
    ("Qin prophage (ydfU)", 1_637_000, 1_638_500),
    ("CPS-53/KpLE1 prophage (yfdM-yfdT)", 2_465_000, 2_470_000),
    ("rrn operon (16S rRNA)", 4_028_500, 4_029_500),
]

# Calls that are reference-genome discrepancies or present in the ancestors.
# Built dynamically below, seeded with the two calls seen in nearly every sample.
FIXED_BACKGROUND = {
    ("NZ_CP009273", 3_446_393, "G"),  # rpsJ synonymous
    ("NZ_CP009273", 4_029_087, "T"),  # 16S rRNA
}

# Chromosomal loci whose disruption is a documented route to colistin or
# beta-lactam resistance; used as a targeted screen.
AMR_LOCI = {
    # colistin / polymyxin
    "pmrA", "pmrB", "basR", "basS", "phoP", "phoQ", "mgrB", "eptA", "eptB",
    "pmrC", "arnA", "arnB", "arnC", "arnD", "arnE", "arnF", "arnT", "lpxA",
    "lpxC", "lpxD", "lpxL", "lpxM", "pagP", "mlaA", "yciM", "lapB",
    # beta-lactam / general efflux and permeability
    "ampC", "ampD", "ampG", "ampR", "acrA", "acrB", "acrR", "acrD", "acrE",
    "acrF", "tolC", "marR", "marA", "marB", "soxR", "soxS", "rob", "envZ",
    "ompR", "ompC", "ompF", "ompA", "mrcA", "mrcB", "ftsI", "pbpA", "dacA",
    "dacB", "mrdA", "lpp", "emrR", "emrA", "emrB", "mdtK", "cpxA", "cpxR",
    "rpoB", "rpoC", "gyrA", "gyrB", "parC", "parE", "nfsA", "nfsB", "ribE",
}

PLASMID_OF_STRAIN = {"PL": "pS2313M_msfGFP", "PM": "pS2313RS_mScarlet", "LM": None}

# The PVB donor carries a pOXA-48 variant lacking the 12.3 kb conjugative
# transfer region downstream of the IS1 element (ssb, mobA, traI and the type IV
# secretion genes), whereas PVI conditions use the intact plasmid. breseq only
# reports this deletion as an unassigned new junction (33164|45514) plus missing
# coverage, because the left boundary sits in the IS1 repeat, so it is
# reconstructed here from the junction evidence, checked against the expected
# plasmid variant, and reported as a DEL call.
PVB_DELETION = {
    "replicon": "NZ_MT441554",
    "junction_left": 33_164,   # last retained base before the deletion
    "junction_right": 45_514,  # first retained base after the deletion
    "gene_name": "[IFA14_RS00230]\u2013[IFA14_RS00325]",
    "gene_product": "pOXA-48 conjugative transfer region: ssb, mobA, traI, "
                    "DotD/TraH lipoprotein, type IV secretion system genes",
}
PVB_JUNCTION_TOLERANCE = 5  # bp
PVB_MIN_FREQUENCY = 0.5     # junction frequency needed to call the variant "deleted"

# Samples dropped before any analysis:
#  - the pilot "old_run" sequencing batch
#  - the mixed KAN+PVI and KAN+PN23 communities, not part of the reported design
#  - three PM colonies picked where a red and a green colony were touching
EXCLUDED_RUNS = {"old_run"}
EXCLUDED_CONDITION_SUBSTRINGS = ("KAN+PVI", "KAN+PN23")
EXCLUDED_SAMPLES = {
    "pOXA48_rep2/M3_PM_clone2",
    "pOXA48_rep2/M5_PM_clone1",
    "pOXA48_rep2/M5_PM_clone2",
}

ANCESTRAL_FREQUENCY_CUT = 0.5


def in_masked(replicon: str, position) -> str | None:
    if replicon != "NZ_CP009273" or pd.isna(position):
        return None
    for name, start, end in MASKED_REGIONS:
        if start <= position <= end:
            return name
    return None


def pvb_deletion_evidence(ev: pd.DataFrame, cov: pd.DataFrame) -> pd.DataFrame:
    """Per-sample evidence for the PVB transfer-region deletion of pOXA-48.

    Returns one row per sample mapped against pOXA-48 with the frequency of the
    diagnostic junction (0 when the junction is absent) and whether breseq also
    reported missing coverage over the region."""
    d = PVB_DELETION
    jc = ev[(ev["kind"] == "JC") & (ev["replicon"] == d["replicon"])].copy()
    jc = jc[jc["side_2_seq_id"].eq(d["replicon"])]
    jc = jc[
        (jc["side_1_position"] - d["junction_left"]).abs().le(PVB_JUNCTION_TOLERANCE)
        & (jc["side_2_position"] - d["junction_right"]).abs().le(PVB_JUNCTION_TOLERANCE)
    ]
    if "reject" in jc:
        jc = jc[jc["reject"].isna()]
    jc["freq"] = pd.to_numeric(jc["frequency"], errors="coerce").fillna(1.0)
    freq = jc.groupby("sample_id")["freq"].max()

    mc = ev[(ev["kind"] == "MC") & (ev["replicon"] == d["replicon"])]
    mc = mc[(mc["start"] <= d["junction_left"] + 1_000) & (mc["end"] >= d["junction_right"] - 100)]
    mc_support = set(mc["sample_id"])

    mapped = sorted(set(cov.loc[cov["replicon_id"] == d["replicon"], "sample_id"]))
    return pd.DataFrame({
        "sample_id": mapped,
        "pOXA48_tra_junction_frequency": [float(freq.get(s, 0.0)) for s in mapped],
        "pOXA48_tra_missing_coverage": [s in mc_support for s in mapped],
    })


def expected_pOXA48_variant(condition) -> str | None:
    cond = str(condition)
    if "PVB" in cond:
        return "deleted"
    if "PVI" in cond:
        return "intact"
    return None


def main() -> None:
    mut = pd.read_csv(OUT / "mutations_raw.tsv", sep="\t", low_memory=False)
    ev = pd.read_csv(OUT / "evidence_raw.tsv", sep="\t", low_memory=False)
    cov = pd.read_csv(OUT / "coverage_raw.tsv", sep="\t")

    mut["frequency"] = pd.to_numeric(mut["frequency"], errors="coerce")

    def keep(df):
        mask = ~df["run"].isin(EXCLUDED_RUNS)
        cond = df["condition"].astype(str)
        for sub in EXCLUDED_CONDITION_SUBSTRINGS:
            mask &= ~cond.str.contains(sub, na=False, regex=False)
        mask &= ~df["sample_id"].isin(EXCLUDED_SAMPLES)
        return df[mask].copy()

    dropped = sorted(set(cov["sample_id"]) - set(keep(cov)["sample_id"]))
    mut, ev, cov = keep(mut), keep(ev), keep(cov)
    print(f"excluded {len(dropped)} samples:")
    for s in dropped:
        print(f"  {s}")

    # ---------------------------------------------------------------- coverage
    meta_cols = ["sample_id", "run", "sample_name", "experiment", "replicate",
                 "vial", "condition", "declared_strain", "sample_type", "timepoint_h"]
    meta = cov[meta_cols].drop_duplicates("sample_id")
    wide = cov.pivot_table(
        index="sample_id", columns="replicon", values="coverage_average", dropna=False
    ).reset_index()
    wide = meta.merge(wide, on="sample_id", how="left")
    for col in ["chromosome", "pOXA48", "PN23", "pS2313M_msfGFP", "pS2313RS_mScarlet"]:
        if col not in wide:
            wide[col] = float("nan")

    wide.columns.name = None
    for col, name in [("pOXA48", "pOXA48_copy_number"), ("PN23", "PN23_copy_number"),
                      ("pS2313M_msfGFP", "msfGFP_copy_number"),
                      ("pS2313RS_mScarlet", "mScarlet_copy_number")]:
        wide[name] = wide[col] / wide["chromosome"]

    # marker plasmid call: > 10x chromosome coverage means the strain carries it
    wide["has_msfGFP"] = wide["msfGFP_copy_number"] > 10
    wide["has_mScarlet"] = wide["mScarlet_copy_number"] > 10
    wide["observed_marker"] = [
        "PL" if g and not s_ else "PM" if s_ and not g else "PL+PM" if g and s_ else "none (LM)"
        for g, s_ in zip(wide["has_msfGFP"], wide["has_mScarlet"])
    ]

    # -------------------------------------------------- auxotrophy verification
    # A knockout shows up either as missing coverage (MC) over the gene, or as a
    # DEL / large_substitution call spanning it, depending on breseq's evidence.
    mc = ev[(ev["kind"] == "MC") & (ev["replicon"] == "NZ_CP009273")].copy()
    mc = mc[["sample_id", "start", "end"]]
    big = mut[
        (mut["replicon"] == "NZ_CP009273")
        & mut["mutation_category"].isin(["large_deletion", "large_substitution"])
    ].copy()
    big["start"] = pd.to_numeric(big["position_start"], errors="coerce")
    big["end"] = pd.to_numeric(big["position_end"], errors="coerce")
    spans = pd.concat([mc, big[["sample_id", "start", "end"]]], ignore_index=True).dropna()

    rows = []
    for sample_id in sorted(set(cov["sample_id"])):
        sub = spans[spans["sample_id"] == sample_id]
        deleted, detail = set(), {}
        for gene, (gstart, gend) in AUX_GENES.items():
            # allow a 60 bp margin: the scar boundaries differ slightly per sample
            hit = sub[(sub["start"] <= gstart + 60) & (sub["end"] >= gend - 60)]
            if len(hit):
                deleted.add(gene)
                detail[gene] = f"{int(hit['start'].min())}-{int(hit['end'].max())}"
        rows.append({"sample_id": sample_id, "deleted_aux_genes": ",".join(sorted(deleted)),
                     "deleted_intervals": ";".join(f"{k}:{v}" for k, v in sorted(detail.items()))})
    aux = pd.DataFrame(rows)
    wide = wide.merge(aux, on="sample_id", how="left")
    wide["deleted_aux_genes"] = wide["deleted_aux_genes"].fillna("")

    def genotype_call(row):
        seen = set(row["deleted_aux_genes"].split(",")) - {""}
        for strain, expected in EXPECTED_KO.items():
            if seen == expected:
                return strain
        return "ambiguous"

    wide["genotype_strain"] = [
        genotype_call({"deleted_aux_genes": g}) for g in wide["deleted_aux_genes"]
    ]

    def marker_matches(strain, marker):
        if strain == "LM":
            return marker == "none (LM)"
        return marker == strain

    # the ancestral clones were mapped against the chromosome only, so their
    # fluorescent-marker status cannot be checked from coverage
    wide["marker_check"] = [
        "not applicable" if run == "ancestrals" or st != "clone" or pd.isna(ds)
        else "match" if marker_matches(ds, om)
        else "MISMATCH"
        for run, st, ds, om in zip(wide["run"], wide["sample_type"],
                                   wide["declared_strain"], wide["observed_marker"])
    ]
    wide["genotype_check"] = [
        "not applicable" if st != "clone" or pd.isna(ds)
        else "match" if gs == ds else "MISMATCH"
        for st, ds, gs in zip(wide["sample_type"], wide["declared_strain"],
                              wide["genotype_strain"])
    ]
    wide["identity_ok"] = [
        (st != "clone") or pd.isna(ds) or (gc == "match" and mk in ("match", "not applicable"))
        for st, ds, gc, mk in zip(
            wide["sample_type"], wide["declared_strain"],
            wide["genotype_check"], wide["marker_check"]
        )
    ]

    # ------------------------------------------ pOXA-48 variant (PVI vs PVB)
    # whole-replicon deletions mark samples in which the plasmid is absent
    absent_pOXA48 = set(mut.loc[
        (mut["kind"] == "DEL") & (mut["position"] == 1)
        & (mut["mutation_category"] == "large_deletion")
        & (mut["replicon"] == PVB_DELETION["replicon"]), "sample_id"])
    pvb = pvb_deletion_evidence(ev, cov)
    wide = wide.merge(pvb, on="sample_id", how="left")

    def tra_region(row):
        if pd.isna(row["pOXA48_tra_junction_frequency"]):
            return "not mapped"
        if row["sample_id"] in absent_pOXA48:
            return "plasmid absent"
        f = row["pOXA48_tra_junction_frequency"]
        if f == 0:
            return "intact"
        return "deleted" if f >= 0.99 else f"deleted ({f:.0%})"

    wide["pOXA48_tra_region"] = wide.apply(tra_region, axis=1)
    wide["pOXA48_expected_variant"] = wide["condition"].map(expected_pOXA48_variant)

    def variant_check(row):
        exp = row["pOXA48_expected_variant"]
        if exp is None or row["pOXA48_tra_region"] in ("not mapped", "plasmid absent"):
            return "not applicable"
        f = row["pOXA48_tra_junction_frequency"]
        observed = "deleted" if f >= PVB_MIN_FREQUENCY else "intact"
        return "match" if observed == exp else "MISMATCH"

    wide["pOXA48_variant_check"] = wide.apply(variant_check, axis=1)
    wide["identity_ok"] &= wide["pOXA48_variant_check"].ne("MISMATCH")

    wide.to_csv(OUT / "sample_identity.tsv", sep="\t", index=False)

    # ------------------------------------------------------ mutation filtering
    # add the PVB transfer-region deletion as an explicit DEL call in every
    # sample whose reads support the diagnostic junction
    d = PVB_DELETION
    del_start, del_end = d["junction_left"] + 1, d["junction_right"] - 1
    carriers_pvb = pvb[pvb["pOXA48_tra_junction_frequency"] > 0]
    sample_meta = cov.drop_duplicates("sample_id").set_index("sample_id", drop=False)
    mode = ev.drop_duplicates("sample_id").set_index("sample_id")["polymorphism_mode"]
    pvb_rows = []
    for sid, f in zip(carriers_pvb["sample_id"], carriers_pvb["pOXA48_tra_junction_frequency"]):
        m = sample_meta.loc[sid]
        pvb_rows.append({
            **{c: m[c] for c in meta_cols}, "polymorphism_mode": bool(mode[sid]),
            "kind": "DEL", "gd_id": "PVB", "replicon": d["replicon"],
            "replicon_name": "pOXA48", "position": del_start,
            "detail": str(del_end - del_start + 1), "frequency": f,
            "mutation_category": "large_deletion",
            "position_start": del_start, "position_end": del_end,
            "gene_name": d["gene_name"], "gene_product": d["gene_product"],
        })
    mut["is_pvb_deletion"] = False
    if pvb_rows:
        extra = pd.DataFrame(pvb_rows)
        extra["is_pvb_deletion"] = True
        mut = pd.concat([mut, extra], ignore_index=True)

    mut["masked_region"] = [in_masked(r, p) for r, p in zip(mut["replicon"], mut["position"])]
    # whole-replicon "deletions" are simply an absent plasmid, not a mutation
    mut["is_absent_replicon"] = (
        (mut["kind"] == "DEL")
        & (mut["position"] == 1)
        & (mut["mutation_category"] == "large_deletion")
    )

    # background = anything fixed in an ancestral clone, plus reference discrepancies
    anc = mut[(mut["run"] == "ancestrals") & (~mut["is_absent_replicon"])].copy()
    anc["freq"] = anc["frequency"].fillna(1.0)
    anc = anc[anc["freq"] >= ANCESTRAL_FREQUENCY_CUT]
    background = set(FIXED_BACKGROUND)
    background |= set(zip(anc["replicon"], anc["position"], anc["detail"]))
    keys = list(zip(mut["replicon"], mut["position"], mut["detail"]))
    mut["is_background"] = [k in background for k in keys]
    # the transfer-region deletion defines the PVB plasmid variant: it is carried
    # by the LM_PVB donor but must still be reported in the evolved samples
    mut.loc[mut["is_pvb_deletion"], "is_background"] = False

    # record which ancestor(s) carry each background call, so that donor-specific
    # variants can be told apart from variants of the plasmid-free parents
    carriers: dict[tuple, list[str]] = {}
    for key, sample in zip(zip(anc["replicon"], anc["position"], anc["detail"]),
                           anc["sample_name"]):
        carriers.setdefault(key, []).append(sample)
    mut["ancestors_carrying"] = [
        ",".join(sorted(set(carriers.get(k, [])))) for k in keys
    ]

    # the engineered auxotrophy knockouts show up as large deletions or
    # substitutions spanning lysA, pheA or metA
    def is_knockout(row):
        if row["replicon"] != "NZ_CP009273":
            return False
        if row["mutation_category"] not in ("large_deletion", "large_substitution"):
            return False
        start = pd.to_numeric(row["position_start"], errors="coerce")
        end = pd.to_numeric(row["position_end"], errors="coerce")
        if pd.isna(start) or pd.isna(end):
            return False
        return any(start <= gs + 60 and end >= ge - 60 for gs, ge in AUX_GENES.values())

    mut["is_knockout"] = mut.apply(is_knockout, axis=1)

    mut["call_class"] = "candidate"
    mut.loc[mut["is_knockout"], "call_class"] = "engineered knockout"
    mut.loc[mut["is_absent_replicon"], "call_class"] = "absent replicon"
    mut.loc[mut["masked_region"].notna() & ~mut["is_knockout"], "call_class"] = "masked repeat/prophage"
    mut.loc[mut["is_background"] & mut["call_class"].eq("candidate"), "call_class"] = "ancestral background"
    # a call seen in an ancestral clone is by definition not a de novo mutation of
    # that clone, whatever its frequency there; the PVB transfer-region deletion
    # is the exception, reported as de novo in the donor as well
    mut.loc[mut["run"].eq("ancestrals") & mut["call_class"].eq("candidate")
            & ~mut["is_pvb_deletion"], "call_class"] = "ancestral background"

    mut.to_csv(OUT / "mutations_classified.tsv", sep="\t", index=False)

    cand = mut[mut["call_class"] == "candidate"].copy()
    cand["amr_locus"] = cand["gene_name"].astype(str).apply(
        lambda g: any(part.strip("[]") in AMR_LOCI for part in str(g).replace("/", " ").split())
    )
    cand.to_csv(OUT / "candidate_mutations.tsv", sep="\t", index=False)

    print("samples:", wide["sample_id"].nunique())
    print("\ncall classes:")
    print(mut["call_class"].value_counts().to_string())
    print("\nidentity check failures:")
    bad = wide[~wide["identity_ok"]]
    print(bad[["sample_id", "declared_strain", "genotype_strain", "deleted_aux_genes",
               "observed_marker", "msfGFP_copy_number", "mScarlet_copy_number"]].round(2).to_string(index=False)
          if len(bad) else "  none")
    print("\npOXA-48 transfer-region deletion (PVB variant) per mapped sample:")
    chk = wide[wide["pOXA48_tra_region"] != "not mapped"]
    print(chk[["sample_id", "condition", "pOXA48_copy_number", "pOXA48_tra_region",
               "pOXA48_tra_missing_coverage", "pOXA48_expected_variant",
               "pOXA48_variant_check"]].round(2).to_string(index=False))

    print("\ncandidate mutations touching an AMR locus:")
    amr = cand[cand["amr_locus"]]
    print(amr[["sample_id", "replicon", "position", "detail", "gene_name",
               "mutation_category", "frequency", "polymorphism_mode"]].to_string(index=False)
          if len(amr) else "  none")


if __name__ == "__main__":
    main()
