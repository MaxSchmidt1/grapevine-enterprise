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
`kirk.py` takes a fastq-file (can be gziped) and checkes every read for the presence of full amplicon sequences of the nine SSR markers and then writes one fastq file and one summary file per maker.
Optionally `--locus` can be used to specify only specific markers to be checked and `--max-mismatches` can be used to specify the number of mismatches per primer (default 2). Furthermore each marker as a specifed size range to which a tolerance window (`--slack`) can be added (default 20). Additionally a cut-off read length can be specified under which reads will be ignored with `--min-read-length` (default 0).


```
kirk.py --input ${WD}/input/${SEED}.fastq.gz --output ${WD}/profiles_new_new/${SEED}/${SEED} --max-mismatches 4 --slack 20 --min-read-length 100
```

## Running `spock.py`
`spock.py` takes the amplicon fastq files and builds a read length matrix for all markers with it while also calling SSR alleles. It has two modes (ambiguous,soft), where ambigous will not attempt to call near homozygous alleles and generally only provide one locus unless the second one is clearly distinct from the first one.

```
$ spock.py -h
usage: spock.py [-h] --basename BASENAME --output OUTPUT [--min-reads MIN_READS] [--stutter-ratio STUTTER_RATIO] [--homozygosity-mode {ambiguous,soft}] [--homo-threshold HOMO_THRESHOLD] [--min-length MIN_LENGTH] [--max-length MAX_LENGTH]
                [--matrix-scale {raw,ce}] [--soft-model-improvement SOFT_MODEL_IMPROVEMENT] [--soft-min-h2-ratio SOFT_MIN_H2_RATIO] [--soft-min-raw-ratio-margin SOFT_MIN_RAW_RATIO_MARGIN] [--soft-cross-parity-markers SOFT_CROSS_PARITY_MARKERS]
                [--soft-bracket-threshold SOFT_BRACKET_THRESHOLD]

Call SSR alleles from extracted amplicon FASTQ files.

options:
  -h, --help            show this help message and exit
  --basename BASENAME   Input basename, e.g. '01' -> reads 01.VVS2.fastq.gz etc.
  --output OUTPUT       Output basename for PNG, TSV, and matrix files.
  --min-reads MIN_READS
                        Minimum total reads per marker to attempt allele calling (default: 1000).
  --stutter-ratio STUTTER_RATIO
                        Peak height fraction below which a +-2/4/6 bp satellite is treated as polymerase stutter and suppressed (default: 0.10). Used in ambiguous mode only; soft mode uses its own stutter model.
  --homozygosity-mode {ambiguous,soft}
                        How to handle close allele pairs (difference <= 4 bp): 'ambiguous' reports only the dominant allele; 'soft' uses a stutter-aware two-peak model to recover the second allele when justified by the data (default: ambiguous).
  --homo-threshold HOMO_THRESHOLD
                        If the second candidate peak is below this fraction of the dominant peak it is treated as homozygous (default: 0.30). Used in ambiguous mode only.
  --min-length MIN_LENGTH
                        Minimum fragment length to include (default: 100).
  --max-length MAX_LENGTH
                        Maximum fragment length to include (default: 320).
  --matrix-scale {raw,ce}
                        Whether the count matrix row indices (sizes) should be on the raw sequencing scale or the CE-calibrated scale (default: raw).
  --soft-model-improvement SOFT_MODEL_IMPROVEMENT
                        [soft mode] Minimum relative drop in residual score that the two-peak model must achieve over the one-peak model (default: 0.25).
  --soft-min-h2-ratio SOFT_MIN_H2_RATIO
                        [soft mode] Minimum fitted height of second allele relative to first allele (default: 0.08).
  --soft-min-raw-ratio-margin SOFT_MIN_RAW_RATIO_MARGIN
                        [soft mode] Required margin by which the candidate's raw ratio must exceed the maximum plausible stutter ratio at that offset (default: 0.05).
  --soft-cross-parity-markers SOFT_CROSS_PARITY_MARKERS
                        [soft mode] Comma-separated list of markers that may carry alleles of either parity (default: VVMD27,VVMD32). Set to an empty string to disable cross-parity candidate search.
  --soft-bracket-threshold SOFT_BRACKET_THRESHOLD
                        [soft mode] If the fitted secondary/dominant height ratio is below this value, the secondary allele is bracketed as low-confidence (default: 0.30).
```

A typical command looks like this

```
spock.py --basename ${WD}/profiles_new_new/${SEED}/${SEED} --output ${WD}/profiles_new_new/${SEED}/${SEED} --homozygosity-mode soft --min-length 100 --max-length 300 --matrix-scale ce
```

## Running `mccoy.py`
`mccoy.py` Uses the SSR profile and compares it against a copy of the VIVC database to assign a variety name if possible.
It requires at least the genotype results from `spock.py` (`--genotype`), the reference database as a tsv-file (`--reference`) and an output file-basename (`--output`).

Optional parameters include the number of displayed hits with `--top-n` (default 10), a mismatch tolerance specifying how much marker-reference sizes can differ from called markers (`--match-tolerance`) (default 1 bp) and the minimum number of matches required to assign a variety call (`--min-anchor-matches`) (default 3).

```
$ mccoy.py 
usage: mccoy.py [-h] --genotype GENOTYPE --reference REFERENCE --output OUTPUT [--top-n TOP_N] [--match-tolerance MATCH_TOLERANCE] [--min-anchor-matches MIN_ANCHOR_MATCHES]
mccoy.py: error: the following arguments are required: --genotype, --reference, --output
[mschmidt@frankfurt ~]$ mccoy.py -h
usage: mccoy.py [-h] --genotype GENOTYPE --reference REFERENCE --output OUTPUT [--top-n TOP_N] [--match-tolerance MATCH_TOLERANCE] [--min-anchor-matches MIN_ANCHOR_MATCHES]

Match a sample SSR genotype against the VIVC reference database.

options:
  -h, --help            show this help message and exit
  --genotype GENOTYPE   Genotype TSV produced by ssr_genotype_sample.py.
  --reference REFERENCE
                        VIVC SSR profiles TSV (vivc_ssr-profiles.tsv).
  --output OUTPUT       Output basename; writes <output>.variety_match.tsv.
  --top-n TOP_N         Number of top candidates to report (default: 10).
  --match-tolerance MATCH_TOLERANCE
                        +-bp tolerance when comparing alleles (default: 1).
  --min-anchor-matches MIN_ANCHOR_MATCHES
                        Minimum number of heterozygous anchor loci that must match (both alleles) for a candidate to be reported (default: 3).
```

A typical command looks like this

```
mccoy.py --genotype ${WD}/profiles_new_new/${SEED}/${SEED}.genotype.tsv --reference ${DATAPATH}/vivc_ssr-profiles.tsv --output ${WD}/profiles_new_new/${SEED}/${SEED}
```
