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
busco -m geno -l $BUSCO/culicidae_odb12 -c 6 -o Anoste_raw -i 
```


