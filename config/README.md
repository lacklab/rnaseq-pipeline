In order to configure your analysis, make changes to `config.yaml`.

# 1. `SAMPLES`
A tab-separated file with the following example should be provided to specify the samples:

| Name             | Condition    |
|------------------|--------------|
| GSM6752476_Vh_1  | Vh           |
| GSM6752477_Vh_2  | Vh           |
| GSM6752478_Vh_3  | Vh           |
| GSM6752479_Dht_1 | Dht          |
| GSM6752480_Dht_2 | Dht          |
| GSM6752481_Dht_3 | Dht          |

Name: Sample name for IP experiment

Control: Input control for ChIPseq

# 2. `UNITS`
A tab-separated file with the following example should be provided to specify the units (all files that need to be preprocessed):

PS: SRR ID units will be fetched from SRA

| Name             | Unit | Fastq1                                                         | Fastq2                                                         | Library     | Dataset     |
|------------------|------|----------------------------------------------------------------|----------------------------------------------------------------|-------------|-------------|
| GSM6752476_Vh_1  | 1    | SRR22381692                                                    | -                                                              | Paired      | GSE218556   |
| GSM6752477_Vh_2  | 1    | SRR22381691                                                    | -                                                              | Paired      | GSE218556   |
| GSM6752478_Vh_3  | 1    | SRR22381690                                                    | -                                                              | Paired      | GSE218556   |
| GSM6752479_Dht_1 | 1    | SRR22381689                                                    | -                                                              | Paired      | GSE218556   |
| GSM6752480_Dht_2 | 1    | SRR22381688                                                    | -                                                              | Paired      | GSE218556   |
| GSM6752481_Dht_3 | 1    | SRR22381687                                                    | -                                                              | Paired      | GSE218556   |

Name: Sample name for processing the data (bam generation)

Unit: Unit no (it will merge the replicates)

Fastq1: Path to the 1st FASTQ file

Fastq2: Path to the 2nd FASTQ file

Dataset: (optional but recommended): The unique dataset identifier.


# 3. `OUTPUT`
- `REF` : The reference genome to be analysed on.
- `RUN`
    - `QC` : Decides if qc analysis will be performed. `True/False` 
    - `QUANT` : Decides if salmon quantification will be performed. This path includes Trimming, STAR alignment and Salmon Quantification. `True/False` 
- `Level` #TODO
    - `Gene`
    - `Transcript`


```
OUTPUT:
    REF: hg19  
    RUN:
        QC: False
        QUANT: TRUE
    LEVEL:
        - Gene
        - Transcript
```

# 4. `REF_hg19` 
- `FA` : Path to genome FASTA file.
- `CHROM_SIZES` : Path to genome's chromosome sizes. 
- `GTF` : Path to gene Annotation file. `GTF` is necessary. 
- `TFA` : Path to transcript FASTA file. (You can create this from previous `GTF` with [gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml).)
- `STAR_IDX` : Path to STAR index. (Generate the STAR index using the files you put here.)

