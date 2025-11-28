# 03_GenomeAssembly

## 00_assemblyRaw

#### Assembling long reads

```
ln -s /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/00_reads/SRR11672506.fastq.gz
conda activate assembly
wtdbg2 -x rs -g 227054799 -t 8 -i SRR11672506.fastq.gz -fo Anoste_raw 
```

**Visualizing assembly output**

```
wtpoa-cns -t 8 -i Anoste_raw.ctg.lay.gz -fo Anoste_raw
```

**Busco Assembly quality check**

```
busco -m geno -l $BUSCO/culicidae_odb12 -c 6 -o Anoste_raw -i Anoste_raw.fasta
```

#### Mapping shot and long reads to start polishing process

**Mapping script**
```
#!/bin/bash
# Mapping short reads on assembly. Mutate SAM into BAM.
# Short reads
minimap2 -ax sr --MD -t 6 Anoste_raw.fasta SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq > Anoste_raw_sr.sam
samtools view -Sb Anoste_raw_sr.sam > Anoste_raw_sr.bam
rm Anoste_raw_sr.sam
samtools sort -@6 -o Anoste_raw_sr_sorted.bam Anoste_raw_sr.bam
samtools index Anoste_raw_sr_sorted.bam
rm Anoste_rae_sr.bam

# Long reads
minimap2 -ax sr --MD -t 6 Anoste_raw.fasta SRR11672506.fastq.gz > Anoste_raw_lr.sam
samtools view -Sb Anoste_raw_lr.sam > Anoste_raw_lr.bam
rm Anoste_raw_lr.sam
samtools sort -@6 -o Anoste_raw_lr_sorted.bam Anoste_raw_lr.bam
samtools index Anoste_raw_lr_sorted.bam
rm Anoste_rae_lr.bam
```
**Run script**
```
conda activate assembly
nano mapping.sh
bash mapping.sh
```
**mosdepth**
```
mosdepth -n --fast-mode --by 500 <...>
```

## 01_polishing

**usign hypo**
```
echo -e "$R1\n$R2" > <READS_PATH>
hypo -d <DRAFT_Contigs> -r @<READS_PATH> -s <APPROXIMATE_GENOMESIZE> -c <SHORT_READSCOVERAGE> -b <SORTED_BAM_SR> -B <SORTED_BAM_PB> -t <NUMBER_THREADS>
```

