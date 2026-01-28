# Genome Annotation and Functional Enrichmet

## Annotation

Functional annotation assigns biological significance to our genomic data. To provide this process there many tools like [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), [diamond](https://github.com/bbuchfink/diamond) and [HIMMER](http://hmmer.org/). 
**BLAST** uses a classic alignment tool; **Diamond** is an accelerated alternative to BLAST, optimized for massive datasets with minimal loss of sensitivity; **HMMER** instead uses Hidden Markov Model (HMM) profiles to detect remote homology. It typically queries profile databases like Pfam (protein families) or the InterPro consortium databases. 

### Input
To functionally explore our data efficienty, we must first generate a clean input file. The optimal strategy is to select the longest protein sequence from each trimmed orthogroup, using `longest_protein_OG.sh` script .
This approach ensures we use the most informative representative for each gene family identified in the orthology analysis. Down below here is reported the correct path for trimmed orthogroups. 

> /00_lab_genomica_comparata/05_OG.Inference_Phylogenomic/04_trimmed/01_prova_trimmed_02_disco_OG$

```bash
bash ../../../99_scripts/longest_protein_OG.sh  
```

### Databases
This a list of databases used in the annotation process executed by **diamond**?

+ **Nr** (Non-redundant): Comprehensive collection from GenPept, Swiss-Prot, PIR, PDB, and RefSeq. Protein content.  
+ **Nt** : Nucleotide sequence collection
+ **Swiss-Prot**: Manually annotated and reviewed proteins (part of UniProt). Protein content 
+ **Pfam**: A large collection of protein families represented by MSAs and HMMs. Domains and Families. 

### Run diamond 

Diamond is optimized for massive datasets, >1 million proteins. Command line provides different flag to manage tool sensitivity settings and resource usage. 

```bash
diamond makedb --in /var/local/diamond_db/nr.gz --db ./nr_diamond
```

-----

## Gene Ontology (GO) Annotation

To link our protein sequences to biological functions, we use Gene Ontology terms. There are several tools available to infer GO terms from protein sequences, such as PANNZER and eggNOG-mapper, both web-based. While effective, these tools can sometimes produce redundant results. For this pipeline, we will use InterProScan, a comprehensive command-line tool that scans sequences against the InterPro member databases including Pfam, PRINTS, SUPERFAMILY, etc.; to ensure robust functional annotation.

```bash
/home/PERSONALE/dbs/interproscan-5.65-97.0/interproscan.sh -i <LONGEST_PROTEINS_INPUT> -goterms -pa -b <OUTPUT-FILE-BASE> -cpu <N_CPUS>
```

> Due to server issues we coudn't be able to run InterProScan, we gave longest_protein file to our Professor who run command line for us. 

-----

## Gene Ontology (GO) Functional enrichment


### Input 

```bash
awk -F'\t' '{
  gsub(/@.*/,"",$1); gsub(/\([^)]*\)/,"",$2); gsub(/\|/,",",$2);
  split($2, a, ",");
  for(i in a) if(a[i]!="") seen[$1,a[i]]=1
}
END {
  for(k in seen){
    split(k, b, SUBSEP)
    groups[b[1]] = (groups[b[1]] ? groups[b[1]] "," b[2] : b[2])
  }
  for(g in groups) print g "\t" groups[g]
}' <(cut -f1,14 longest_pietro.tsv) | grep -v "-" > go_back.ts
```
