# grapevine-enterprise
This is a pipeline consisting of three scripts (`kirk.py`, `spock.py`, `mccoy.py`) to identify grapevine varieties from Oxford Nanopore Sequencing amplicons of the nine standard grapevine SSR-markers (VVS2, VVMD5, VVMD7, VVMD25, VVMD27, VVMD28, VVMD32, VrZAG62, VrZAG79)

`kirk.py` will scan the input fastq files and extract the amplicon sequences including the primer sequences themselves and create one binned fastq-file for each marker.
`spock.py` will use the fastq file to generate a size ditribution count matrix of all reads between 100 and 300 bp for all markers in a samples and try to call the peaks in this size distribution and produce a genotype.tsv file with the identified markers
`mccoy.py` will compare this against a tsv database file containing all varieties annotated in the VIVC database (as of April 2026) and produce a variety_match.tsv file highlighting the most likely variety identification

## Data Preparation
Best Run life basecalling and demultiplexing during sequencing so only the individual fastqs per barcode need to be concatanated e.g.
```
for i in {01..96}
  do
    cat /mnt/raw/nanopore/20260326_Previti_SSR_02/Run01/20260326_1220_P2S_00089-B_PBE58034_3737082c/fastq_pass/barcode${i}/*.fastq.gz input/${i}.fastq.gz
  done
```
(I used super high-accuracy basecalling and set the Q-score cut-off at Q20, which removes more than half of the reads but we do not need a lot of data)

## Running `kirk.py`


kirk.py --input ${WD}/input/${SEED}.fastq.gz --output ${WD}/profiles_new_new/${SEED}/${SEED} --max-mismatches 4 --slack 20 --min-read-length 100
