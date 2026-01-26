# Orthology Inference 

First step to begin a comparative genomic pipeline, after creating our [dataset](https://github.com/bacconipietro/00_lab_genomica_comparata/blob/main/00_practice/00_data/README.md), is making a phylogenetic inference on dataset species. We must identify **Orthologs**: genes separated by speciation and distinguish them from **Paralogs**: genes separated by duplication. Then we need to order ONLY ortholog genes in **Orthogroups**, sets of genes descended from a single gene in the last common ancestor of a group of species. All these complex inferences based on analyzing thousands of genes are simply reachable running an effective tool as **Orthofinder**.  
It uses a gene-centric graph clustering approach: 

+ Identify all orthologs and paralogs.
+ Infer unrooted gene trees for every orthogroup.
+ Reconstruct a robust species tree from the data.

## Orthofinder

**Running orthofinder**
```bash
ln -s /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/03_dataset/02_proteome/00_aa/Anoste.faa . #This is the assembly we made in precedent directories
conda activate orthofinder
#[orthofinder]
orthofinder -t 8 -a 8 -f 02_raw_proteoms/ #Giving in input all annotated genomes
```
In the end software gives in output a folders network with a serie of different results, which will be used to process phylogenetic inference. Run statistics are reported here `Statistics_Overall.tsv`. 

-----

## Paralogy filtering - DISCO

As we sad before we need to separate and delete paralogs from orthologs. We use a tool called [DISCO](https://github.com/jsdoublel/DISCO/blob/master/README.md), it decomposes complex groups into smaller, distinct sub-clusters that are strictly consistent with orthology and it requires a specific FASTA header format to correctly identify species. Fortunately, the syntax we established earlier. 

```bash
conda activate tree
cd Orthofinder/Results_Dec01/Resolved_Gene_Trees
#[tree]
while IFS=' ' read -r OG tree; do python3 ~/Lab_genomica_comparata/00_practice/99_scripts/disco.py -i <(echo "$tree") -o ../../../01_DISCO/${OG/:/}.nwk -d "|" -m 4 --remove_in_paralogs --keep-labels --verbose >> ../../../01_DISCO/disco.log; done < <(sed -E 's/[A-Z][a-z]{5}_//g; s/\)n[0-9]*+/\)/g' Resolved_Gene_Trees.txt)
```
There's a chance to find empty files, we delete them.  
```bash
find . -size 0 -print > empty_disco.txt
find . -size 0 -delete
```

After DISCO has decomposed the complex orthogroups, we need to generate individual sequence files for each new sub-cluster. This was achieved using the `split_disco_output.sh` script, which correlates the phylogenetic tree partitions .nwk with the corresponding amino acid data. The script effectively splits the bulk dataset into separate protein files .faa ready for multiple sequence alignment.

Here is the script:
```bash
#!/bin/bash

echo "Splitting DISCO outputs ..."

mkdir split

# Split the multiline output in multiple outputs
for disco in *.nwk; do
	OGname=$(basename -s .nwk "$disco")
	split -d -a 2 -l 1 --additional-suffix=.nwk "$disco" split/"$OGname"_
done

mkdir ../disco_OG

# Define original sequence folder with $1
OG_folder=$(realpath $1)

# recreate orthogroups using original ones. It is MANDATORY to substitute all '-' with '_' because it is te only way to grap exactly
# the meant sequence.

echo "Recreating orthogroups ..."

cd split

for tree in *_*.nwk; do
    name="${tree%.nwk}"
    OG="${name%%_*}"  # Extract the part before the number
    output="../../disco_OG/$name.faa"

    # Extract unique sequence names from tree
    grep -o -E "[A-Z][a-z]{5}\.[^:]+" "$tree" | sort -u > tmp_seqs.txt

    # Extract sequences in one go using awk
    awk -v seqs=tmp_seqs.txt '
        BEGIN {
            while ((getline < seqs) > 0) wanted[$1] = 1
        }
        /^>/ {
            seq = substr($0, 2)
            keep = (seq in wanted)
        }
        keep { print }
    ' "$OG_folder/$OG.faa" >> "$output"
done

# Clean up temporary file
rm -f tmp_seqs.txt
```

**Run:**
```bash
cd 05_OG.Inference_Phylogenomic/01_DISCO
bash ../../99_scripts/split_disco_output.sh ../OrthoFinder/Results_Dec01/Orthogroup_Sequences #In input we need to give extant orthogroups sequences 
mv disco_OG 02_disco_OG
```

-----

## Alignments and trimming 
Following the decomposition of orthogroups, the pipeline proceeds to Multiple Sequence Alignment (MSA) and Trimming. While DISCO identifies which sequences belong together, MSA determines how they line up. This establishes the column-by-column positional information necessary to distinguish evolutionary signal from noise. We need to align SIngle Copy Orthologue Sequences, from this folder we exctract 200 random sequence to proceed.  

```bash
mkdir 05_OG.Inference_Phylogenomic/03_aligned 05_OG.Inference_Phylogenomic/04_trimmed
cd 05_OG.Inference_Phylogenomic/Orthofinder/Results_Dec01/Single_Copy_Orthologue_Sequences
ls *.fa | shuf -n 200 > species_tree_OG.txt
```

We run [mafft](https://github.com/GSLBiotech/mafft) for **MSA**. 
```bash
conda activate sequence
#[sequence]
for OG in $(cat species_tree_OG.txt); do mafft --auto --anysymbol "$OG" > ../../../03_aligned/${OG/.fa/_aligned.fa}; done
```

We run [BMGE](https://github.com/PurdueRCAC/Biocontainers/blob/main/docs/source/bmge/bmge.rst) for **Trimming**. It calculates the entropy (variability) at each position and removes sites that are too variable (saturated) to provide reliable evolutionary information, discarding regions with excessive gaps or ambiguous alignments. 
```bash
cd 03_aligned
mkdir 00_single_complete
mv *_aligned.fa 00_single_complete
conda deactivate    #you need to deactivate all the open environments to run 'bmge'
for OG in *; do bmge -i "$OG" -t AA -m BLOSUM62 -e 0.5 -g 0.4 -of ../../04_trimmed/${0G/_aligned.fa/_trimmed.fa}; done
```
-----

## Species tree Inference 

In the final stage it's possible to execute phylogenetic analysis. At first we need to **concat** all trimmed alignments using [AMAS](https://github.com/marekborowiec/AMAS/blob/master/amas/AMAS.py). Before runnig the concatenation we must delete noise patterns and making sure the script works on same headers. 
 
```bash
cd 04_trimmed
sed -i.old -E 's/\|.+$//' *   # sed to delete undesire patterns
../../../99_scripts/AMAS.py concat -y nexus -i *.fa -f fasta -d aa -t conc_species_tree   #concatenation
```

After concatenation we are able to run [IQtree](https://iqtree.github.io/doc/). 
```bash
mkdir 05_OG.Inference_Phylogenomic/05_tree/00_species
conda activate tree
#[tree]
iqtree -m TESTNEW -b 100 -s conc_species_tree --prefix 05_tree/00_species/species_tree -nt 9
```
