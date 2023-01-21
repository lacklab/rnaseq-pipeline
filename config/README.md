In order to configure your analysis, make changes to `config.yaml`.

# 1. `SAMPLES`
A tab-separated file with the following example should be provided to specify the samples:

| Name            | Conrol       |
|-----------------|--------------|
| AR_NKI_1_N      | -            |
| AR_NKI_7_N      | -            |
| LuCAP176nlAR    | LuCAP176nlIN |
| LuCAP176nlFOXA1 | LuCAP176nlIN |
| LuCAP176wzAR    | LuCAP176wzIN |
| LuCAP176wzFOXA1 | LuCAP176wzIN |

Name: Sample name for IP experiment

Control: Input control for ChIPseq

# 2. `UNITS`
A tab-separated file with the following example should be provided to specify the units (all files that need to be preprocessed):

PS: SRR ID units will be fetched from SRA

| Name             | Unit | Fastq1                                                         | Fastq2                                                         | Library     |
|------------------|------|----------------------------------------------------------------|----------------------------------------------------------------|-------------|
| AR_NKI_1_N       | 1    | SRR11856198                                                    | -                                                              | Single      |
| AR_NKI_7_N       | 1    | SRR11856199                                                    | -                                                              | Single      |
| LuCAP176nlAR     | 1    | raw-data/LuCAP176nlAR_CKDL220029757-1A_HLYT2DSX5_L2_1.fq.gz    | raw-data/LuCAP176nlAR_CKDL220029757-1A_HLYT2DSX5_L2_2.fq.gz    | Paired      |
| LuCAP176nlFOXA1  | 1    | raw-data/LuCAP176nlFOXA1_CKDL220029753-1A_HLYT2DSX5_L2_1.fq.gz | raw-data/LuCAP176nlFOXA1_CKDL220029753-1A_HLYT2DSX5_L2_1.fq.gz | Paired      |
| LuCAP176nlIN     | 1    | raw-data/LuCAP176nlIN_CKDL220029757-1A_HLYT2DSX5_L2_1.fq.gz    | raw-data/LuCAP176nlIN_CKDL220029757-1A_HLYT2DSX5_L2_2.fq.gz    | Paired      |
| LuCAP176wzAR     | 1    | raw-data/LuCAP176wzAR_CKDL220029755-1A_HLYT2DSX5_L2_1.fq.gz    | raw-data/LuCAP176wzAR_CKDL220029755-1A_HLYT2DSX5_L2_2.fq.gz    | Paired      |
| LuCAP176wzFOXA1  | 1    | raw-data/LuCAP176wzFOXA1_CKDL220029756-1A_HLYT2DSX5_L2_1.fq.gz | raw-data/LuCAP176wzFOXA1_CKDL220029756-1A_HLYT2DSX5_L2_2.fq.gz | Paired      |   
| LuCAP176wzIN     | 1    | raw-data/LuCAP176wzIN_CKDL220029756-1A_HLYT2DSX5_L2_1.fq.gz    | raw-data/LuCAP176wzIN_CKDL220029756-1A_HLYT2DSX5_L2_2.fq.gz    | Paired      |

Name: Sample name for processing the data (bam generation)

Unit: Unit no (it will merge the replicates)

Fastq1: Path to the 1st FASTQ file

Fastq2: Path to the 2nd FASTQ file


# 3. `OUTPUT`
- `RUN`
    - `QC` : Decides if qc analysis will be performed. `True/False` 
    - `PEAKS` : Decides if peak calling will be performed. `True/False` 
    - `BWS` : Decides if bigwig files will be generated. `True/False` 
- `BW_NORMALIZATIONS` : (If `BWS:True`) normalizes the bigwig files accordingly. `rawcounts/RPM`
- `BAMPROCESS_PARAMS` : `samtools view` parameters to filter read accordingly. PS: This pipeline is uniqe reads and paired samples aware, thus `-f` and `-F` used builtin, but you can change MAPQ threshold.
- `MACS_THRESHOLD` : The significancy treshold for `macs3 callpeak -q`.

```
OUTPUT:
    RUN:
        QC: False
        PEAKS: True
        BWS: True
    BW_NORMALIZATIONS:
        - rawcount
    BAMPROCESS_PARAMS: -q 30
    MACS_THRESHOLD: 0.01 
```

# 4. `REF`
   - `Name` : Name of the genome
   - `FA`: Path to the FASTA file of the reference genome
   - `BWA_IDX`: Path to the BWA index files (prefix)
   - `CHROM_SIZES` : Path to the chrom.sizes file of the reference genome

  