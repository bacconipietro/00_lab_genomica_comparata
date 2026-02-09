# Genome Annotation and Functional Enrichmet

## Annotation

Functional annotation assigns biological significance to our genomic data. To provide this process there many tools like [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), [diamond](https://github.com/bbuchfink/diamond) and [HIMMER](http://hmmer.org/). BLAST uses a classic alignment tool; Diamond is an accelerated alternative to BLAST, optimized for massive datasets with minimal loss of sensitivity; HMMER instead uses Hidden Markov Model (HMM) profiles to detect remote homology. It typically queries profile databases like Pfam (protein families) or the InterPro consortium databases. 

For this pipeline we chose **Double Indexing Alignment of Next-generation Sequencing Data** (DIAMOND) which acts as the high-speed engine that compares sequences against massive reference databases, maintaining high level of sensitivity.  

### Input
To functionally explore our data efficienty, we must first generate a clean input file. The optimal strategy is to select the longest protein sequence from each trimmed orthogroup, using `longest_protein_OG.sh` script .
This approach ensures we use the most informative representative for each gene family identified in the orthology analysis. Down below here is reported the correct path for trimmed orthogroups. 

> /00_lab_genomica_comparata/05_OG.Inference_Phylogenomic/04_trimmed/trimmed_02_disco_OG

```bash
cd /00_lab_genomica_comparata/05_OG.Inference_Phylogenomic/04_trimmed/trimmed_02_disco_OG
bash ../../../99_scripts/longest_protein_OG.sh
less longest_protein_OGs.txt
```

### Databases
This a list of databases used in the annotation process executed by **Diamond**?

+ **Nr** (Non-redundant): Comprehensive collection from GenPept, Swiss-Prot, PIR, PDB, and RefSeq. Protein content.  
+ **Nt** : Nucleotide sequence collection
+ **Swiss-Prot**: Manually annotated and reviewed proteins (part of UniProt). Protein content 
+ **Pfam**: A large collection of protein families represented by MSAs and HMMs. Domains and Families. 

## Diamond 

Diamond is optimized for massive datasets, >1 million proteins and it's used to find homologous sequences for our Orthogroups, which serves as the foundation for assigning functional metadata like protein names. Command line provides different flag to manage tool sensitivity settings and resource usage. 

```bash
diamond makedb --in /home/STUDENTI/pietro.bacconi/00_lab_genomica_comparata/05_OG.Inference_Phylogenomic/04_trimmed/trimmed_02_disco_OG/longest_protein_OGs.txt --db ./nr_diamond
```

> Due to server issues we coudn't be able to run diamond, we gave longest_protein file to our Professor who run command line for us.

-----

## Gene Ontology (GO) Annotation

To link our protein sequences to biological functions, we use Gene Ontology terms. There are several tools available to infer GO terms from protein sequences, such as PANNZER and eggNOG-mapper, both web-based. While effective, these tools can sometimes produce redundant results. For this pipeline, we will use InterProScan, a comprehensive command-line tool that scans sequences against the InterPro member databases including Pfam, PRINTS, SUPERFAMILY, etc.; to ensure robust functional annotation.

```bash
/home/PERSONALE/dbs/interproscan-5.65-97.0/interproscan.sh -i longest_protein_OGs.txt -goterms -pa -b longest_pietro.tsv -cpu <N_CPUS>
```

> Due to server issues we coudn't be able to run InterProScan, we gave longest_protein file to our Professor who run command line for us. 

-----

## Gene Ontology (GO) Functional enrichment

The next step is to understand the biological meaning of these changes. We employ **Gene Ontology (GO) Enrichment Analysis**. This is a robust statistical framework designed to translate lists of gene identifiers into recognizable biological themes. The core objective is to identify functional "signals" that rise above the background "noise" of the entire genome. We compare a **background** dataset, composed by every orthogroup that has been assigned at least one GO term, with a **foreground** one, the subset of "interesting" genes, identified by CAFE as significantly expanded, which are suspected to be under positive selection or involved in a specific phenotypic trait. 


#### Input

+ At first we need to remove extra metadata, pipe separators, and descriptions from `longest_pietro.tsv` which are not allowed in R script.  
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
}' <(cut -f1,14 longest_pietro.tsv) | grep -v "-" > go_back.tsv
```

+ Before running the enrichment, you need a list of orthogroups you want to test. In this workflow, these are extracted from the CAFE results and put in `interesting.txt`. The grep is wrote based on biological constraints that are relevant for our analysis. To investigate the genomic basis of cold adaptation, we specifically targeted gene families showing significant evolutionary turnover in two key lineages: Toxrut and Anosin. These species were selected because they exhibit phenotypic traits associated with cold tolerance or inhabit temperate ranges, unlike their tropical relatives.

```bash
grep "Toxrut<1>\*" Gamma_asr.tre | grep "<3>\*" | grep "<2>_" | grep "<4>_" | grep "<5>_" | grep "<6>_" > CAFE_relevant_2_OG.txt

nano interesting.txt

OG0000027
OG0000045
OG0000060
OG0000102
OG0000121
```
> This passage must be run in the CAFE directory that we selected as best model previously. In this case the correct path is:
00_lab_genomica_comparata/07_GeneFamilies_Evolution/00_1L/4K/9N/Gamma_asr.tre 

-----

#### R script 
Here is reported R script to run the enrichment. Results are avaible in `/p2_Enrichment` directory. 
```R
library(tidyverse)
library(topGO)

gene_universe <- readMappings(file =
                                "go_back_collapsed.tsv)
geneUniverse <- names(gene_universe)

genesOfInterest <- read.table("interesting.txt",header=FALSE)
list_interest1 <- list( "name_interest" = genesOfInterest)

#upload of gene of interest
GOenrichment <- function(trait, trait_name) {
  if (!dir.exists("01_enrichment")) {
    dir.create("01_enrichment")
  }
  
  genesOfInterest <- as.character(trait[[1]]) #as vector not character 
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  
  print(trait_name)
  
  ontology_values = c("BP", "MF", "CC")
  
  GOdata_list <- lapply(ontology_values, function(ontology_value) {
    GOdata_name <- paste("GOdata_", ontology_value, sep = "")
    # annot = annFUN.gene2GO this imparts the program which annotation it should use. In this case, it is specified that it will be in gene2GO format and provided by the user.
    # gene2GO = gene_universe is the argument used to tell where is the annotation
    assign(GOdata_name, new("topGOdata", ontology=ontology_value, allGenes=geneList, annot = annFUN.gene2GO, gene2GO = gene_universe))
  })
  
  elim_list <- lapply(seq_along(ontology_values), function(i) {
    elim_name <- paste("elim_", ontology_values[i], sep="")
    assign(elim_name, runTest(GOdata_list[[i]], algorithm="elim", statistic="fisher"))
  })
  
  results_elim <- function(GO_data, elim_data) {
    num_nodes <- min(1000, length(elim_data@score))
    resulte <- GenTable(GO_data, Classic_Fisher = elim_data,
                        orderBy = "Classic_Fisher", topNodes=num_nodes, numChar=1000)
    resulte$Classic_Fisher <- as.numeric(resulte$Classic_Fisher)
    resulte <- subset(resulte, Classic_Fisher < 0.05)
    return(resulte)
  }
  
  results_elim_list <- lapply(seq_along(ontology_values), function(i) {
    resulte_name <- paste("resulte_", ontology_values[i], sep="")
    assign(resulte_name, envir = .GlobalEnv, results_elim(GOdata_list[[i]], elim_list[[i]]))
  })
  
  write_elim_results <- function(result, ontology_value, trait_name) {
    table_name <- paste("01_enrichment/topGOe_", trait_name, "_", ontology_value, ".txt", sep="")
    write.table(result, file=table_name, quote=F, sep = "\t", row.names = F)
  }
  
  lapply(seq_along(ontology_values), function(i) {
    write_elim_results(results_elim_list[[i]], ontology_values[i], trait_name)
  })
}

GMT <- function(trait, trait_name) {
  genesOfInterest <- as.character(trait$V1) #as vector not character 
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  
  print(paste("Working on", trait_name, sep = " "))
  
  ontology_values = c("BP", "MF", "CC")
  
  GOdata_list <- lapply(ontology_values, function(ontology_value) {
    GOdata_name <- paste("GOdata_", ontology_value, sep = "")
    # annot = annFUN.gene2GO this imparts the program which annotation it should use. In this case, it is specified that it will be in gene2GO format and provided by the user.
    # gene2GO = gene_universe is the argument used to tell where is the annotation
    assign(GOdata_name, new("topGOdata", ontology=ontology_value, allGenes=geneList, annot = annFUN.gene2GO, gene2GO = gene_universe))
  })
  
  elim_list <- lapply(seq_along(ontology_values), function(i) {
    elim_name <- paste("elim_", ontology_values[i], sep="")
    assign(elim_name, runTest(GOdata_list[[i]], algorithm="elim", statistic="fisher"))
  })
  
  results_elim <- function(GO_data, elim_data) {
    resulte <- GenTable(GO_data, Classic_Fisher = elim_data, orderBy = "Classic_Fisher", topNodes=1000, numChar=1000)
    resulte$Classic_Fisher <- as.numeric(resulte$Classic_Fisher)
    resulte <- subset(resulte, Classic_Fisher < 0.05)
    return(resulte)
  }
  
  results_elim_list <- lapply(seq_along(ontology_values), function(i) {
    resulte_name <- paste("resulte_", ontology_values[i], sep="")
    assign(resulte_name, envir = .GlobalEnv, results_elim(GOdata_list[[i]], elim_list[[i]]))
  })
  
  #transform in GMT format
  list_OG_GO <- function(GO_term, GOdata, gene_of_interest){
    genes <- intersect(genesInTerm(GOdata, GO_term)[[1]], gene_of_interest)
    return(paste(genes, collapse = ","))
  }
  
  result_GTM_list <- lapply(seq_along(ontology_values), function(i) {
    results_elim_list[[i]]$Genes <- unlist(sapply(results_elim_list[[i]]$GO.ID, function(GO_term){list_OG_GO(GO_term, GOdata_list[[i]], genesOfInterest)}))
    return(results_elim_list[[i]])
  })
  
  print("Done! writing results.")
  
  write_gtm_results <- function(result, ontology_value, trait_name) {
    table_name <- paste("02_enrichment/GTM_", trait_name, "_", ontology_value, ".gtm", sep="")
    write.table(result, file=table_name, quote=F, sep = "\t", row.names = F)
  }
  
  lapply(seq_along(ontology_values), function(i) {
    write_gtm_results(result_GTM_list[[i]], ontology_values[i], trait_name)
  })
}

#Complete function to perform enrichment
##this particular syntax has been necessary since it was impossible to give the function the trait name it was computing.
GO_enrichment <- function(list) {
  lapply(seq_along(list), function(i) {
  GOenrichment(list[[i]], names(list)[i])
  })
}

GMT_parsing <- function(list) {
  lapply(seq_along(list), function(i) {
    GMT(list[[i]], names(list)[i])
  })
}

#Run the complete function
GO_enrichment(social_not_aculeata)
GMT_parsing(list_notext_11)

#search for specific gene
genesOfInterest <- as.character(notext_biggest$V1) #as vector not character 
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
topGO <- new("topGOdata", ontology="BP", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = gene_universe)
print(intersect(genesInTerm(topGO, "GO:0051289")$'GO:0051289', notext_biggest$V1)) #to investigate which genes are annotated with a particular GO term
```
