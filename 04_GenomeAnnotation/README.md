# Genome Annotation

This phase focuses on identifying the physical locations of genomic features within the assembly. The result of this process is usually a .gff (General Feature Format) or .gtf file, text files that contain the coordinates of every feature relative to your assembly. Genome functional annotation yields  descriptions of specific biological functions for each nucleotide sequence. 

-----

## Repeats elements and transposons masking

It's necessary to exclude repetitive regions and transposons from gene annotation because genome prediction software could trace false positive signals, mistaking trasposons protein for host genes. With this purpose we **mask** these elements. To identify repetitive elements, we utilize a standard workflow involving [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler/blob/master/README.md) and [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker/blob/master/README.md). 
+ **RepeatModeler** scans the genome de novo to build a library of consensus repeat sequences specific to this organism. 
+ **RepeatMasker** uses that library to search the genome and "mask" (hide) the repeats so they don't interfere with gene prediction.

Since RepeatModeler is computationally intensive, we have already received the consensus library. It's avaible here:

> /home/PERSONALE/mirko.martini3/01_2024/00_Data/02_annotation/Anoste_RepeatModeler_library.fa

We will run RepeatMasker internally as a step within the MAKER pipeline and summarize the results at the end. 

-----

##  First Round - MAKER

We use MAKER for the genome annotation. Unlike tools that take a long string of command-line flags, MAKER relies on specific control files to manage its inputs and parameters. 
This command generates three specific configuration files required for the pipeline. MAKER relies on these files instead of long command-line arguments.
```bash
maker -CTL
```

+ `maker_exe.ctl` contains the paths to all the executable tools (BLAST, RepeatMasker, etc.). MAKER usually finds these automatically.

+ `maker_bopts.ctl` defines filtering thresholds for BLAST and other alignment algorithms. We will maintain the default parameters for this file.

+ `maker_opts.ctl` is the primary configuration file containing input data paths and run settings. This is the only file we need to modify for our analysis.


```bash
    #-----Genome (these are always required)
    genome= <GENOME> #genome sequence (fasta file or fasta embeded in GFF3 file)
    organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic
    
    #-----Re-annotation Using MAKER Derived GFF3
    maker_gff= #MAKER derived GFF3 file
    est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
    altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
    protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
    rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
    model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
    pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
    other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no
    
    #-----EST Evidence (for best results provide a file for at least one)
    est= #set of ESTs or assembled mRNA-seq in fasta format
    altest= #EST/cDNA sequence file in fasta format from an alternate organism
    est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
    altest_gff= #aligned ESTs from a closly relate species in GFF3 format
    
    #-----Protein Homology Evidence (for best results provide a file for at least one)
    protein= <PROTEOME> #protein sequence file in fasta format (i.e. from mutiple oransisms)
    protein_gff=  #aligned protein homology evidence from an external GFF3 file
    
    #-----Repeat Masking (leave values blank to skip repeat masking)
    model_org=<EMPTY> #select a model organism for RepBase masking in RepeatMasker
    rmlib= <RepeatModeler_library> #provide an organism specific repeat library in fasta format for RepeatMasker
    repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
    rm_gff= #pre-identified repeat elements from an external GFF3 file
    prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
    softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)
    
    #-----Gene Prediction
    snaphmm= #SNAP HMM file
    gmhmm= #GeneMark HMM file
    augustus_species= #Augustus gene prediction species model
    fgenesh_par_file= #FGENESH parameter file
    pred_gff= #ab-initio predictions from an external GFF3 file
    model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
    est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
    protein2genome=<0 OR 1> #infer predictions from protein homology, 1 = yes, 0 = no
    trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
    snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
    unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
    
    #-----Other Annotation Feature Types (features MAKER doesn't recognize)
    other_gff= #extra features to pass-through to final MAKER generated GFF3 file
    
    #-----External Application Behavior Options
    alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
    cpus=<CPUS> #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)
    
    #-----MAKER Behavior Options
    max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
    min_contig=1 #skip genome contigs below this length (under 10kb are often useless)
    
    pred_flank=200 #flank for extending evidence clusters sent to gene predictors
    pred_stats=<0 or 1> #report AED and QI statistics for all predictions as well as models
    AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
    min_protein=50 #require at least this many amino acids in predicted proteins
    alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
    always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
    map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
    keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)
    
    split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
    single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
    single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
    correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes
    
    tries=2 #number of times to try a contig if there is a failure for some reason
    clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
    clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
    TMP= #specify a directory other than the system default temporary directory for temporary files
```
Commit this changes and run the code line.

+ `genome=:` specify the path to your genome assembly FASTA file.
+ `protein=:` Provide the path to external proteome evidence. 
+ `model_org=:` Remove the default value (leave this blank to avoid biasing the repeat masking).
+ `rmlib=:` Enter the path to the custom RepeatModeler library generated in the previous step.
+ `protein2genome=1:` Change from 0 to 1. This enables MAKER to infer gene predictions directly from protein homology (essential since we are providing external proteome evidence).
+ `cpus=:` Set threads number.  
+ `phred_stats=1:` Change from 0 to 1 to enable the calculation of Phred quality scores for each annotation.
+ `min_protein=50:` Set the minimum length in amino acids for a predicted protein to be included in the final annotation.
+ `alt_splice=1:` Change from 0 to 1 to allow the prediction of alternative splicing isoforms.
+ `split_hit=:` Remove default value 

```bash
maker -base <OUTPUT PREFIX>
```

Proceed merging the resulting .fasta and .gff files.
```bash
fasta_merge -d <DATASTORE INDEX FILE>
gff3_merge -d <DATASTORE INDEX FILE>
```
Conclude the first round with statistic evaluation using [AGAT](https://github.com/NBISweden/AGAT). Additionally, to verify the first round proterome it's possible run BUSCO.     
```bash
#[GAAS]
agat_sp_statistics.pl --gff file.gff -o <output_file>
agat_sq_repeats_analyzer.pl -i <input_file> -o <output_file>
```

-----

## Second Round - SNAP and Augustus 

### SNAP 
**SNAP** is a probabilistic gene predictor designed to identify protein-coding genes, especially in novel genome assemblies. It utilizes Hidden Markov Models (HMMs) to recognize gene features like start/stop codons and splice sites. This predictor is highly effective when trained on organism-specific data. To achieve this training, we rely on Fathom (included in the SNAP suite). Fathom processes the initial gene evidence (from MAKER or other sources) to create a high-quality training set, which is then used to teach SNAP how to recognize genes in our specific genome.

```bash
#I suggest you to create a small script with all these commands because are really fast
maker2zff -c 0 -e 0  -l 80 -x 0.1 -d <datastore_index> #To extract gene models based on mutiple filter criterion. It transforms maker output into .zff files, the format for SNAP
fathom <.ANN FILE> <.DNA FILE> -gene-stats #Print some summary statistics of the selected gene models
fathom <.ANN FILE> <.DNA FILE> -validate #Validate gene models and print summary statistics. output is in the stdout and cannot be redirected
fathom <.ANN FILE> <.DNA FILE> -categorize 1000 #Extract gene models together with 1000 bp at both ends for training
fathom <UNI.ANN FILE> <UNI.DNA FILE> -export 1000 -plus #Export and convert unique genes to strand plus. Tipically, only the unique genes are sued for training and testing. The value limits the intergenic sequence at the ends.
forge <export.ann> <export.dna> #call the parameter estimation program, better in another directory
hmm-assembler.pl <NAME> <FORGE DIRECTORY> > <OUTPUT HMM FILE>
```

### Augustus
**Augustus** is another accurate gene predictor available for eukaryotic genomes. It uses a Generalized Hidden Markov Model (HMM) to identify gene structures such as exons, introns, UTRs without requiring external transcript data. he optimal way to train Augustus is to curate a strict dataset of high-confidence gene models. This involves using tools like AGAT, but it could be computationally intensive. For this reason we chose an easier path using BUSCO, which automatically finds conserved genes and uses them to train the Augustus. 

```bash
#[sequence]
busco -i ../../03_GenomeAssembly/03_scaffolding/Anoste_chr.fasta -c 30 -l $BUSCO/culicidae_odb12 --augustus --long -m genome --out Anoste_cu --augustus_parameters='--progress=true'
```

### Run second round annotation after training 

Create a copy of your previous configuration file and modify the following parameters:

+ `protein_gff / rm_gff:` Add the paths to the alignment GFF files generated in the previous round, contains protein alignment and repeat masking evidence.
+ `protein:` Clear this field as we are now using the GFFs above. 
+ `snaphmm:` Add the path to your trained SNAP HMM file.
+ `augustus_species:` Set this to your specific species name (e.g., Anoste_cu).
+ `protein2genome / est2genome`: Set both to 0 (Disabled). We are now relying on the trained ab initio predictors rather than direct evidence alignment.
+ `pred_stats:` Ensure this remains set to 1 (Enabled) to calculate quality scores, especially if planning a third annotation round.

Execute MAKER again and collect the resulting GFF, protein, and transcript files, it's possible to carry out multiple rounds. It's reasonable to record a significant increase in annotation quality when moving from round 1 to round 2. 

### Evaluation of gene set  

This command line returns round annotation summary statistics.
```bash
#[GAAS]
gaas_maker_merge_outputs_from_datastore.pl Anoste_rnd2.maker.output/ #script for summary statistic, merge fastas and gffs
agat_sq_repeats_analyzer.pl --gff maker_mix.gff -o Anoste_rnd2_repeats.txt #use maker_mix, which is the complete one
```

Best advice to verify quality annotation is to compare summary results with bibliography knowledge, especially mean gene length, mean exon length, mean number of exon and introns per gene.
