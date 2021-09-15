# ICGC-ARGO Open-access Regions
ICGC-ARGO Open-access regions definitions and bed files

## Definitions for Open-access Regions 
- Genomic elements are defined according to the definitions in the table based on [GENCODE annotations v38](https://www.gencodegenes.org/human/release_38.html)
- Open-access = CDS + Exon + UTR + Promoter + splice site + small ncRNA
- Control-access = complementary regions of Open-access upon whole genome


| Genomic Elements           | Definition                                                                                                                                                                                                                                                                     | Operation                                                                                                                                                                                                                                                                                                        |
|---------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Coding elements (CDS)      | The set of coding bases collapsed across all coding transcripts with a given GENCODE gene ID.                                                                                                                                                                                  | feature_type=CDS                                                                                                                                                                                                                                                                                                 |
| Exon                       | The set of exon bases collapsed across all transcripts with a given GENCODE gene ID.                                                                                                                                                                                           | feature_type=exon                                                                                                                                                                                                                                                                                                |
| 5'UTR (UTR5)               | The set of 5’UTR bases collapsed across all coding transcripts with a given GENCODE gene ID.                                                                                                                                                                                   | feature_type=UTR                                                                                                                                                                                                                                                                                                 |
| 3'UTR (UTR3)               | The set of 3’UTR bases collapsed across all coding transcripts with a given GENCODE gene ID.                                                                                                                                                                                   | feature_type=UTR                                                                                                                                                                                                                                                                                                 |
| Protein-coding promoter    | Regions extending 200 bases in upstreaming direction from all protein-coding transcripts’ transcription start sites (5’ ends). <br />Bases were collapsed across all coding transcripts with a given GENCODE gene ID.                                                                | feature_type=transcript<br />transcript_type=protein_coding<br />locate the TSS for each transcript<br />get the regions (-200, 0) upstreaming of TSS                                                                                                                                                                               |
| Protein-coding splice site | Intronic regions extending 6 bases from donor splice sites and 20 bases from acceptor splice sites were collected for all coding transcripts<br />Bases were collapsed across all coding transcripts with a given GENCODE gene ID.                                                 | feature_type=transcript<br />  transcript_type=protein_coding<br />locate all the introns and donor/acceptor splice sites for each transcript<br />get the regions (-3,6) and (-20,3) around donor/acceptor splice sites |
| LncRNA                     | lncRNA transcripts were defined based on annotations from GENCODE. <br />The elements were made by collapsing bases across transcripts with given GENCODE gene ID.                                                                                                                   | feature_type=transcript<br /> transcript_type=lncRNA                                                                                                                                                                                                                                                                  |
| LncRNA promoter            | Regions extending 200 bases in upstreaming direction from all lncRNA transcripts’ transcription start sites (5’ ends). <br />Bases were collapsed across all lncRNA transcripts with a given gene GENCODE ID.                                                                        | feature_type=transcript<br /> transcript_type=lncRNA<br />locate the TSS for each transcript<br />get the regions (-200,0) upstreaming of TSS                                                                                                                                                                                        |
| LncRNA splice site         | Intronic regions extending six bases from donor splice sites and 20 bases from acceptor splice sites were collected for all lncRNA transcripts. <br />Bases were collapsed across all lncRNA transcripts with a given gene GENCODE ID.                                               | feature_type=transcript<br /> transcript_type=lncRNA<br />locate all the introns and donor/acceptor splice sites for each transcript<br />get the regions (-3, 6) and (-20, 3) around donor/acceptor splice sites          |
| smallRNA                   | smallRNA transcripts were defined based on GENCODE transcripts with biotype: ['Mt_rRNA', 'Mt_tRNA', 'miRNA', 'misc_RNA', 'rRNA', 'scRNA', 'snRNA', 'snoRNA', 'ribozyme', 'sRNA', 'scaRNA']. <br />Bases were collapsed across all smallRNA transcripts with a given gene GENCODE ID. | feature_type=transcript<br />transcript_type in ['Mt_rRNA', 'Mt_tRNA', 'miRNA', 'misc_RNA', 'rRNA', 'scRNA', 'snRNA', 'snoRNA', 'ribozyme', 'sRNA', 'scaRNA']                                                                                                                                                        |
| smallRNA promoter          | Regions extending 200 bases in upstreaming direction from all smallRNA transcripts’ transcription start sites (5’ ends). <br />Bases were collapsed across all smallRNA transcripts with a given gene GENCODE ID.                                                                    | feature_type=transcript<br />transcript_type in ['Mt_rRNA', 'Mt_tRNA', 'miRNA', 'misc_RNA', 'rRNA', 'scRNA', 'snRNA', 'snoRNA', 'ribozyme', 'sRNA', 'scaRNA']<br />locate the TSS for each transcript<br />get the regions (-200, 0) upstreaming of TSS                                                                              |


## Script for generating Open-access bed files
### Use default setting
```
cd script
./gtf_region.py
```

### Script usage
```
./gtf_region.py -h
usage: gtf_region.py [-h] [-g GTF_DIR] [-u GTF_URL] [-v GTF_VERSION] [-l GLEN]
                     [-w WINDOW_SIZE] [-b BED_DIR] [-r REGIONS [REGIONS ...]]
                     [-o OPEN_ACCESS_REGIONS [OPEN_ACCESS_REGIONS ...]]

optional arguments:
  -h, --help            show this help message and exit
  -g GTF_DIR, --gtf_dir GTF_DIR
                        local directory for GTF file: only ENSEMBL or GENCODE GTF file accepted
  -u GTF_URL, --gtf_url GTF_URL
                        specify URL to download GTF file. Use GENCODE by default
  -v GTF_VERSION, --gtf_version GTF_VERSION
                        specify GTF release version to download
  -l GLEN, --glen GLEN  file for genome size
  -w WINDOW_SIZE, --window WINDOW_SIZE
                        the value of w in calculating TSS regions. Default: 200
  -b BED_DIR, --bed_dir BED_DIR
                        directory for output bed files
  -r REGIONS [REGIONS ...], --regions REGIONS [REGIONS ...]
                        specify regions to generate bed files
  -o OPEN_ACCESS_REGIONS [OPEN_ACCESS_REGIONS ...], --open_access_regions OPEN_ACCESS_REGIONS [OPEN_ACCESS_REGIONS ...]
                        specify regions to be included into open access tier
```



## Bed files for Open-access Regions
- You can find the generated bed files for Open-access Regions at your local folder `data/hg38/bed/gencode.v**/`. E.g.,
```
data/hg38/bed
└── gencode.v38
    ├── cds.gencode_v38.20210915.bed.gz
    ├── control_access.gencode_v38.20210915.bed.gz
    ├── exon.gencode_v38.20210915.bed.gz
    ├── intron.gencode_v38.20210915.bed.gz
    ├── lncRNA.gencode_v38.20210915.bed.gz
    ├── lncRNA_promoter.gencode_v38.20210915.bed.gz
    ├── lncRNA_splice_site.gencode_v38.20210915.bed.gz
    ├── open_access.gencode_v38.20210915.bed.gz
    ├── protein_coding.gencode_v38.20210915.bed.gz
    ├── protein_coding_promoter.gencode_v38.20210915.bed.gz
    ├── protein_coding_splice_site.gencode_v38.20210915.bed.gz
    ├── smallRNA.gencode_v38.20210915.bed.gz
    ├── smallRNA_promoter.gencode_v38.20210915.bed.gz
    ├── smallRNA_splice_site.gencode_v38.20210915.bed.gz
    ├── utr3.gencode_v38.20210915.bed.gz
    └── utr5.gencode_v38.20210915.bed.gz
```

- Here each genomic element has its own bed file for example: `exon.gencode_v**.{date}.bed.gz`, `smallRNA.gencode_v**.{date}.bed.gz`, `utr3.gencode_v**.{date}.bed.gz` and etc.
- `open_access.gencode_v**.{date}.bed.gz` contains the union of all the Open-access regions defined above.
- `control_access.gencode_v**.{date}.bed.gz` contains the complementary regions of Open-access upon whole genome.
- However ONLY `open_access.gencode_v**.{date}.bed.gz` are kept in git.




