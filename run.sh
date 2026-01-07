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

# run breseq, using 12 cores concurrently
parallel --jobs 12 --progress < jobs.sh
