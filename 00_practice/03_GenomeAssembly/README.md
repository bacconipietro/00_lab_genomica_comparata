# 03_GenomeAssembly

## 00_assemblyRaw

### Assembling long reads

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

### Mapping short and long reads to start polishing process

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
**Calculate short reads coverage with mosdepth** 
```
conda activate assembly
mosdepth -n --fast-mode --by 500 Anoste_raw_sr Anoste_raw_sr_sorted.bam
zcat Anoste_raw_sr.regions.bed.gz | awk '{sum += $4;count++} END {print sum / count}' > <OUTFILE_coverage>
```

## 01_polishing

**Using hypo to polish raw assembly**
```
echo -e "$R1\n$R2" > Sr.path
conda activate assemply
hypo -d Anoste_raw.fasta -r @Sr.path -s 227054799 -c 136 -b Anoste_raw_sr.bam -B Anoste_raw_lr.bam -t 6
mv hypo_Anoste_raw.fasta Anoste_pol.fasta 
```
### Quality check on polished assembly**
  **N50**
```
assembly-stats Anoste_pol.fasta
```
  **Busco**
```
fold -w 80 Anoste_pol.fasta Anoste_pol_one.fasta
conda activate sequence
busco -m geno -l $BUSCO/culicidae_odb12 -c 8 -o Anoste_pol_busco -i Anoste_pol.fasta
cat short_summary #reading busco result summery
```
  **KAT**
```
conda activate kat
kat comp -t 8 lo Anoste_pol 'SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq' Anoste_pol.fasta
```
## 02_contaminants

### Necessary steps before assembly decontamination
**Re-mapping short reads**
```
minimap2 --secondary=no --MD -ax sr -t 8 Anoste_pol.fasta SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq | samtools view -Sb - > Anoste_pol_sr.bam
samtools sort -@8 -o Anoste_pol_sr_sorted.bam Anoste_pol_sr.bam
rm Anoste_pol_sr.bam
samtools index Anoste_pol_sr_sorted.bam
```
**Taxonomic contig annotation**
```
conda activate assembly
blastn -query <ASSEMBLY> -db <PATH/TO/nt/> -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -max_target_seqs 25 -max_hsps 1 -num_threads 25 -evalue 1e-25 -out <OUTFILE>  #we didn't run this lines

####mancano dei nomi
```
### Decontamination and assembly visualitation
```
blobtools create -i Anoste_pol.fasta -b <MAPPING_FILE> -t <BLASTN_FILE> -o <OUTPUT_PREFIX>
blobtools view -i <JSON_FILE> -o Anoste
blobtools plot -i <JSON_FILE> -o Anoste

####da finire mancano dei nomi
```
**Different command to visualize result**
```
grep -v '^##' Anoste_blobDB_table.txt | column -t -s $'\t' | le ss`  ##example
```
**Grep and Awk ONLY Arthropoda contings**
```
grep "Arthropoda" Anoste.Anoste_blob.blobDB.table.txt > contig_arthropoda.tsv
wc -l contig_arthropode.tsv
grep -w -A1 "ctg1" Anoste_pol.fasta | head
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' Anoste_pol.fasta | grep -w  -Ff <(cut -f1 contig_arthropoda.tsv) - | tr "\t" "\n" > Anoste_decontaminated.fasta
```

### Scaffolding 
```
ragtag.py correct -t 20 <REFERENCE_GENOME> <DRAFT_GENOME>
ragtag.py scaffold -C -t 20 -o <OUTPUT_DIX> <REFERENCE_GENOME> <CORRECTED_DRAFTGENOME>
```
