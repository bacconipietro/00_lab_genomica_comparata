# Comparative Genomics

Welcome to Comparative Genomics github repository. Here are reported all directories that resumes an enitre workflow course of Unibo, University of Bologna, Master degree in Biodiversity and Evolution. 
During this semester we focused on learning assembly and comparative inference skills using lots of bioinformatic tools. This repository documents a dual-phase genomic analysis centered on the invasive malaria vector *Anopheles stephensi*. We conducted a genome assembly on *A. stephensi* and a multi-species study involving 5 distinct Culicidae genomes selected from diverse thermal niches. By contrasting temperate lineages against tropical counterparts, we aim to detect genomic signatures of cold adaptation. 
Every folder reports a specific passage, click on titles to figure out what we have done!  


+ [00_data](https://github.com/bacconipietro/00_lab_genomica_comparata/tree/main/00_practice/00_data): 
  This section details the construction of a comparative genomic dataset comprising *Anopheles stephensi* and five other mosquito species to investigate cold adaptation, starting with automated NCBI data retrieval based on     climatic niches. It documents the subsequent curation pipeline, including filtering for longest isoforms using AGAT, removing pseudogenes with internal stop codons, and standardizing sequence headers for downstream         orthology analysis.

+ [01_ComputerEnvs](https://github.com/bacconipietro/00_lab_genomica_comparata/tree/main/01_ComputerEnv):
  Defines the computational infrastructure, detailing the specific Conda environments configured to isolate dependencies for genome assembly, sequence processing, and phylogenetic analysis.  

+ [02_KmerBased_GenomeSurvey](https://github.com/bacconipietro/00_lab_genomica_comparata/tree/main/02_KmerBased_GenomeSurvey):
  Outlines the pre-assembly quality control workflow, detailing the use of FastQC for initial assessment and Trimmomatic to remove low-quality bases and adapters from raw FASTQ reads. It describes the subsequent genome survey using K-mer analysis (KAT) and GenomeScope to estimate genomic characteristics like heterozygosity and size prior to assembly.

+ [03_GenomeAssembly](https://github.com/bacconipietro/00_lab_genomica_comparata/tree/main/03_GenomeAssembly):
  This section details the hierarchical assembly pipeline followed by iterative quality assessments using N50, BUSCO, and KAT to ensure structural and biological integrity. It outlines the subsequent refinement stages, including error polishing, taxonomic decontamination to filter non-Arthropoda sequences, and final reference-based scaffolding. 

+ [04_GenomeAnnotation](https://github.com/bacconipietro/00_lab_genomica_comparata/tree/main/04_GenomeAnnotation):
  This section outlines the genome annotation workflow followed by an initial evidence-based MAKER round to map protein homology and generate training data. Here is explained the iterative refinement process, which involves training ab-initio gene predictors on high-confidence models and executing a second MAKER pass to produce the final structural annotation, validated via AGAT statistics.

+ [05_OG.Inference_Phylogenomic](https://github.com/bacconipietro/00_lab_genomica_comparata/tree/main/05_OG.Inference_Phylogenomic):
  This section describes the phylogenetic workflow, starting with the inference of orthologs and paralogs using OrthoFinder to cluster genes into Orthogroups, followed by the decomposition of complex multi-copy families into strictly orthologous sub-clusters using DISCO. It details the subsequent species tree reconstruction, which involves selecting single-copy orthologs, performing multiple sequence alignment and trimming, and finally inferring the phylogeny using a supermatrix approach with IQ-TREE.

+ [06_DivergenceTime_Estimation](https://github.com/bacconipietro/00_lab_genomica_comparata/tree/main/06_DivergenceTime_Estimation):
  It explains the construction of an ultrametric time tree using IQ-TREE with the Least Square Dating (LSD2) algorithm, chosen as a faster alternative to Bayesian methods for dating phylogenomic inferences. Calibration process is broke down, using divergence dates retrieved from [TimeTree.org](https://timetree.org/) to constrain specific nodes and the final execution of the molecular clock analysis to estimate divergence times across the phylogeny.

+ [07_GeneFamilies_Evolution](https://github.com/bacconipietro/00_lab_genomica_comparata/tree/main/06_DivergenceTime_Estimation):
  The section details the analysis of gene family turnover using CAFE, describing the preparation of input data and the execution of stochastic birth-death models to estimate global ($1\lambda$) and local ($2\lambda$) evolutionary rates. Here is listed the statistical model selection pipeline, utilizing custom Bash scripts to calculate AIC and BIC scores from likelihood outputs to identify the optimal parameter configuration for detecting significant expansions and contractions.

+ [09_GeneAnnotation_functional_enrichment](https://github.com/bacconipietro/00_lab_genomica_comparata/tree/main/09_GeneAnnotation_functional_enrichment):
  This section describes the functional annotation pipeline, which assigns biological metadata to orthogroups by aligning the longest representative proteins against comprehensive databases. It details the subsequent Gene Ontology (GO) enrichment analysis, employing a custom R script to identify statistically over-represented biological functions in gene families showing significant evolutionary turnover.

-----

## Github features
```
git add <input_file> #[upload a file in the staging area]
git commit -m "" #create a commit
git push #push files from staging area to remote repository
git status #check files which have not been uploaded yet
git mv <any_file> #move any files
git rm <any_file> #remove any files
```

## Screen

```
tmux new -s <screen_name> #create a new screen
tmux a -t <screen_name> #reopen extant screens
tmux ls #list all existing screens
Ctrl+a+d #close screen
Ctrl+d #kill screen
Ctrl+a+c #create a new window in the current session
Ctrl+a+n #move to the window on the right (next) in the current session
Ctrl+a+p #move to the window on the left (previous) in the current session
Ctrl+a+| #split vertically the current window
Ctrl+a+_ #split horizontally the current window
Ctrl+a+[ #enter in read mode
```


