# 02_KmerBased_GenomeSurvey

### Quality check

```bash
# [assembly]
fastqc fastqc SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz
```

### Trimming process

```
conda activate assembly
trimmomatic PE -threads 20 -phred33 SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz SRR11672503_1_paired.fastq SRR11672503_1_unpaired.fastq SRR11672503_2_paired.fastq SRR11672503_2_unpaired.fastq ILLUMINACLIP:/opt/miniforge3/envs/assembly/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> stats_trimmomatic.log
```


### Kat, analizing k-mer distributions

```
conda activate kat
kat hist -t 11 -m 227054799 -o kat_stats SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq
```

