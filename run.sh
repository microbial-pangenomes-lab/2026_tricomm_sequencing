# use the provided conda env
# mamba env create -f environment.yml
# conda activate tricomm

# Download BW25113 genomes
ncbi-genome-download -H -F gff,fasta,genbank,protein-fasta -p 36 -o reference bacteria -t 67989
# get the one we care about
zcat reference/human_readable/refseq/bacteria/Escherichia/coli/K-12_substr._BW25113/GCF_000750555.1_ASM75055v1_genomic.gbff.gz > BW25113.gbk

# prepare out directories
mkdir -p out
mkdir -p out/ancestrals
mkdir -p out/pOXA48_rep1/
mkdir -p out/pOXA48_rep2/
mkdir -p out/PN23_rep1/
mkdir -p out/PN23_rep2/

# run breseq, using 12 cores concurrently
parallel --jobs 12 --progress < jobs.sh

# per-base coverage of the ancestral clones (used by coverage_ancestrals.ipynb)
mosdepth -t 1 out/ancestrals/PL1/coverage out/ancestrals/PL1/data/reference.bam
mosdepth -t 1 out/ancestrals/PL2/coverage out/ancestrals/PL2/data/reference.bam
mosdepth -t 1 out/ancestrals/PM1/coverage out/ancestrals/PM1/data/reference.bam
mosdepth -t 1 out/ancestrals/PM2/coverage out/ancestrals/PM2/data/reference.bam
mosdepth -t 1 out/ancestrals/LM1/coverage out/ancestrals/LM1/data/reference.bam
mosdepth -t 1 out/ancestrals/LM2/coverage out/ancestrals/LM2/data/reference.bam

# parse the breseq output, classify mutations and verify strain identity,
python parse_gd.py
python classify.py
python make_tables.py
