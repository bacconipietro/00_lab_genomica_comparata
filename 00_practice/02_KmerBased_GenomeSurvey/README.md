# Kmer Based Genome Survey

## Quality check for reads

FASTQ is standard text-file format that stores both a nucleotide sequence and its corresponding quality scores. 
We used **_fastqc_** tool to obtain reads quality. Results are reported in two .html files.  


```
# [assembly]
conda activate assembly
fastqc SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz
```


## Trimming 

Trimming is a data cleaning process performed on raw FASTQ files before analysis. It involves computationally removing unwanted sequences from the ends of the reads to ensure that only high-quality data is used for alignment or assembly. We used **_trimmomatic-0.40_** in _assembly_ enviroment. 

```
# [assembly]
trimmomatic PE -threads 8 -phred33 SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz SRR11672503_1_paired.fastq SRR11672503_1_unpaired.fastq SRR11672503_2_paired.fastq SRR11672503_2_unpaired.fastq ILLUMINACLIP:/opt/miniforge3/envs/assembly/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> stats_trimmomatic.log
```

Trimming outputs in this case are: two files for paired reads, two files for unpaired and one .log file for trimmomatic stats. 
It would be a good idea to redo _fastqc_ after trimming process, just to be sure about the quality.   


## K-mer frequency

KAT is a specialized suite designed to evaluate the structure and completeness of genome assemblies. By comparing k-mer distributions between raw sequencing reads and the final assembly, KAT helps verify that the assembly accurately represents the original data. A k-mer is simply a substring of DNA of length _k_. We used **_kat_** command in _kat_ enviroment. 


```
#[kat]
conda activate kat
kat hist -t 8 -m 27 -o Anoste SRR11672503_1_paired_fastqc.html SRR11672503_2_paired_fastqc.html
```

Two main outputs: **Anoste.hist** and **Anoste.png**. To visualize our main results it's necessary to remove all comments and upload Anoste.hist on [genomescope](http://genomescope.org/genomescope2.0/). With this tool we can read coverage results about k-mer distributions, image below. Usually the tallest pick is for homozygosity while the short one rapresents eterozigosity. (Unfortunately our genoscope wasn't correct, for this reason i didn't upload .jpg files). 

