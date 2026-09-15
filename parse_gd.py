#!/usr/bin/env python3
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

SEQ = Path(__file__).resolve().parent
OUT = Path(__file__).resolve().parent

MUT_TYPES = {"SNP", "SUB", "DEL", "INS", "MOB", "AMP", "CON", "INV"}
EVIDENCE_TYPES = {"RA", "MC", "JC", "UN"}

REPLICONS = {
    "NZ_CP009273": "chromosome",
    "NZ_MT441554": "pOXA48",
    "MG557852": "PN23",
    "pS2313M_msfGFP": "pS2313M_msfGFP",
    "pS2313RS_mScarle": "pS2313RS_mScarlet",
    "pS2313RS_mScarlet": "pS2313RS_mScarlet",
}

# BW25113 coordinates of the three auxotrophy genes (from the repo notebooks)
AUX_GENES = {
    "lysA": (2_970_996, 2_972_258),
    "pheA": (2_731_104, 2_732_264),
    "metA": (4_204_208, 4_205_137),
}

# Expected double knockout per TriComm member
EXPECTED_KO = {
    "PL": {"pheA", "lysA"},
    "PM": {"pheA", "metA"},
    "LM": {"lysA", "metA"},
}

POXA48_MAP = {
    "M0": "TriComm KAN no drug",
    "M1": "TriComm KAN 16xMIC AMP",
    "M2": "TriComm KAN+PVI no drug",
    "M3": "TriComm KAN+PVI 16xMIC AMP",
    "M4": "TriComm PVI no drug",
    "M5": "TriComm PVI 16xMIC AMP",
    "M6": "TriComm PVB no drug",
    "M7": "TriComm PVB 16xMIC AMP",
}
PN23_MAP_REP1 = {
    "M1": "TriComm KAN+PN23 2xMIC COL",
    "M2": "TriComm KAN no drug",
    "M3": "TriComm KAN 1xMIC COL",
    "M4": "TriComm KAN 2xMIC COL",
    "M5": "TriComm PN23 no drug",
    "M6": "TriComm PN23 1xMIC COL",
    "M7": "TriComm PN23 2xMIC COL",
}
PN23_MAP_REP2 = {
    "M0": "TriComm KAN no drug",
    "M1": "TriComm PN23 no drug",
    "M2": "TriComm KAN 1xMIC COL",
    "M3": "TriComm PN23 1xMIC COL",
    "M4": "TriComm KAN 2xMIC COL",
    "M5": "TriComm PN23 2xMIC COL",
}


def parse_gd(path: Path) -> tuple[list[dict], list[dict]]:
    """Return (mutations, evidence) rows from one annotated.gd file."""
    mutations, evidence = [], []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            kind = parts[0]
            attrs = {}
            for field in parts[4:]:
                if "=" in field:
                    key, value = field.split("=", 1)
                    attrs[key] = value
            if kind in MUT_TYPES:
                row = {
                    "kind": kind,
                    "gd_id": parts[1],
                    "replicon": parts[3],
                    "position": int(parts[4]) if parts[4].isdigit() else None,
                    "detail": parts[5] if len(parts) > 5 else "",
                }
                row.update(attrs)
                mutations.append(row)
            elif kind in EVIDENCE_TYPES:
                row = {"kind": kind, "gd_id": parts[1], "replicon": parts[3]}
                if kind in {"MC", "UN"}:
                    try:
                        row["start"] = int(parts[4])
                        row["end"] = int(parts[5])
                    except (ValueError, IndexError):
                        continue
                elif kind == "JC":
                    # new junction: (seq_id, position, strand) for both sides
                    try:
                        row["side_1_position"] = int(parts[4])
                        row["side_1_strand"] = int(parts[5])
                        row["side_2_seq_id"] = parts[6]
                        row["side_2_position"] = int(parts[7])
                        row["side_2_strand"] = int(parts[8])
                    except (ValueError, IndexError):
                        continue
                row.update(attrs)
                evidence.append(row)
    return mutations, evidence


def sample_metadata(run: str, sample: str) -> dict:
    """Derive experiment / replicate / sample-type metadata from folder names."""
    meta = {"run": run, "sample_name": sample}
    if run == "ancestrals":
        meta["experiment"] = "ancestral"
        meta["vial"] = None
        meta["declared_strain"] = sample[:2]
        meta["sample_type"] = "clone"
        meta["timepoint_h"] = 0.0
        if "_" in sample:  # LM_PVB: the plasmid-carrying donor clone
            meta["replicate"] = None
            meta["condition"] = f"ancestral clone, {sample.split('_', 1)[1]} donor"
        else:
            meta["replicate"] = sample[2:]
            meta["condition"] = "ancestral clone"
        return meta
    if run == "old_run":
        meta["experiment"] = "pilot"
        meta["replicate"] = None
        meta["condition"] = sample
        meta["vial"] = None
        meta["declared_strain"] = sample.split("_")[0]
        meta["sample_type"] = "population"
        m = re.search(r"_(\d+)h$", sample)
        meta["timepoint_h"] = float(m.group(1)) if m else None
        return meta

    experiment, replicate = run.split("_")
    meta["experiment"] = experiment
    meta["replicate"] = replicate
    vial = sample.split("_")[0]
    meta["vial"] = vial
    if experiment == "PN23":
        table = PN23_MAP_REP1 if replicate == "rep1" else PN23_MAP_REP2
    else:
        table = POXA48_MAP
    meta["condition"] = table.get(vial)
    strain = None
    for candidate in ("PL", "PM", "LM"):
        if f"_{candidate}_" in sample or sample.endswith(f"_{candidate}"):
            strain = candidate
    meta["declared_strain"] = strain
    meta["sample_type"] = "clone" if "clone" in sample else "population"
    m = re.search(r"_(\d+)h", sample)
    meta["timepoint_h"] = float(m.group(1)) if m else None
    return meta


def main() -> None:
    mut_rows, ev_rows, cov_rows = [], [], []
    for summary in sorted(SEQ.glob("out/*/*/data/summary.json")):
        sample_dir = summary.parent.parent
        run = sample_dir.parent.name
        sample = sample_dir.name
        meta = sample_metadata(run, sample)
        meta["sample_id"] = f"{run}/{sample}"

        with summary.open() as fh:
            data = json.load(fh)
        refs = data.get("references", {}).get("reference", {})
        for ref_name, ref in refs.items():
            cov_rows.append(
                {
                    **meta,
                    "replicon_id": ref_name,
                    "replicon": REPLICONS.get(ref_name, ref_name),
                    "length": ref.get("length"),
                    "coverage_average": ref.get("coverage_average"),
                    "reads_mapped": ref.get("num_reads_mapped_to_reference"),
                }
            )
        gd = sample_dir / "data" / "annotated.gd"
        muts, evs = parse_gd(gd)
        # breseq polymorphism mode (-p) is recorded in the command line
        polymorphism = " -p " in gd.read_text(errors="ignore").split("\n")[3]
        for row in muts:
            mut_rows.append({**meta, "polymorphism_mode": polymorphism, **row})
        for row in evs:
            ev_rows.append({**meta, "polymorphism_mode": polymorphism, **row})

    mutations = pd.DataFrame(mut_rows)
    evidence = pd.DataFrame(ev_rows)
    coverage = pd.DataFrame(cov_rows)

    mutations["replicon_name"] = mutations["replicon"].map(REPLICONS).fillna(
        mutations["replicon"]
    )
    evidence["replicon_name"] = evidence["replicon"].map(REPLICONS).fillna(
        evidence["replicon"]
    )

    mutations.to_csv(OUT / "mutations_raw.tsv", sep="\t", index=False)
    evidence.to_csv(OUT / "evidence_raw.tsv", sep="\t", index=False)
    coverage.to_csv(OUT / "coverage_raw.tsv", sep="\t", index=False)
    print(f"samples: {coverage['sample_id'].nunique()}")
    print(f"mutation calls: {len(mutations)}")
    print(f"evidence rows: {len(evidence)}")
    print(mutations["kind"].value_counts().to_string())


if __name__ == "__main__":
    main()
