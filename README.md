# 2026_tricomm_sequencing

Variant calling (breseq) and downstream analysis of TriComm evolution experiments.

## Contents

- `*.gbk`: reference sequences (BW25113 chromosome is downloaded by `run.sh`; plasmids pOXA48, PN23, pS2313M_msfGFP, pS2313RS_mScarlet are provided).
- `out/<experiment>/<sample>/data/`: breseq outputs (`annotated.gd`, `summary.json`) for each sample. Experiments: `ancestrals`, `pOXA48_rep1`, `pOXA48_rep2`, `PN23_rep1`, `PN23_rep2`, `old_run`.
- `out/ancestrals/<clone>/coverage.per-base.bed.gz`: mosdepth per-base coverage of the ancestral clones.
- `run.sh`, `jobs.sh`, `old_run_jobs.sh`: the commands used to run breseq and mosdepth.
- `parse_gd.py`, `classify.py`, `make_tables.py`: parse breseq outputs, classify mutations, build supplementary tables.
- `collate_references.py`: per-replicon coverage summary (`references_collated.csv`).
- `visualize_breseq.py`, `run_visualization_batch.py`: per-sample mutation maps.
- `*.ipynb`: figures.

## Environment

```
mamba env create -f environment.yml
conda activate tricomm
```

## Reproducing the analysis

The breseq and mosdepth steps (`run.sh` / `jobs.sh`) are included for reference only: the read
paths point to our workstation, and their outputs are already in `out/`.
Start from those:

```
python parse_gd.py        # -> mutations_raw.tsv, evidence_raw.tsv, coverage_raw.tsv
python classify.py        # -> mutations_classified.tsv, candidate_mutations.tsv, sample_identity.tsv
python make_tables.py     # -> supplementary_sequencing_tables.xlsx, table_S*.tsv
python collate_references.py  # -> references_collated.csv
python run_visualization_batch.py  # -> plots/
```

Then run the notebooks:

- `coverage.ipynb`: uses `references_collated.csv`.
- `PVB_PVI.ipynb`: uses `sample_identity.tsv`.
- `coverage_ancestrals.ipynb`: uses `out/ancestrals/*/coverage.per-base.bed.gz`.
