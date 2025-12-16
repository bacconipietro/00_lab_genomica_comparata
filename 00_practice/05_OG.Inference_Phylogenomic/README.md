# Orthofinder Results
```
mv Orthofinder ../../05_OG.Inference_Phylogenomic
mkdir 05_OG.Inference_Phylogenomic/01_DISCO
```

## Paralogy filtering 
```
conda activate tree
cd Orthofinder/Results_Dec01/Resolved_Gene_Trees
while IFS=' ' read -r OG tree; do python3 ~/Lab_genomica_comparata/00_practice/99_scripts/disco.py -i <(echo "$tree") -o ../../../01_DISCO/${OG/:/}.nwk -d "|" -m 4 --remove_in_paralogs --keep-labels --verbose >> ../../../01_DISCO/disco.log; done < <(sed -E 's/[A-Z][a-z]{5}_//g; s/\)n[0-9]*+/\)/g' Resolved_Gene_Trees.txt)
```

## Delete empty files
```
find . -size 0 -print > empty_disco.txt
find . -size 0 -delete
```

## Run split_disco script
```
cd 05_OG.Inference_Phylogenomic/01_DISCO
bash ../../99_scripts/split_disco_output.sh ../OrthoFinder/Results_Dec01/Orthogroup_Sequences
mv disco_OG 02_disco_OG
```

## Alignments and trimming Single Copy Orthologue Sequences
```
mkdir 05_OG.Inference_Phylogenomic/03_aligned 05_OG.Inference_Phylogenomic/04_trimmed
cd 05_OG.Inference_Phylogenomic/Orthofinder/Results_Dec01/Single_Copy_Orthologue_Sequences
ls *.fa | shuf -n 200 > species_tree_OG.txt
```

Alignment
```
conda activate sequence
for OG in $(cat species_tree_OG.txt); do mafft --auto --anysymbol "$OG" > ../../../03_aligned/${OG/.fa/_aligned.fa}; done
```

Trimming
```
cd 03_aligned
mkdir 00_single_complete
mv *_aligned.fa 00_single_complete
conda deactivate    #you need to deactivate all the open environments to run 'bmge'
for OG in *; do bmge -i "$OG" -t AA -m BLOSUM62 -e 0.5 -g 0.4 -of ../../04_trimmed/${0G/_aligned.fa/_trimmed.fa}; done
```

## Species tree Inference 

Modify headers
```
cd 04_trimmed
mkdir 00_single_complete
mv *_trimmed.fa 00_single_complete
cd 00_single_complete
sed -i.old -E 's/\|.+$//' *
```

Concat alignments
```
../../../99_scripts/AMAS.py concat -y nexus -i *.fa -f fasta -d aa -t conc_species_tree
```

Run iqtree
```
mkdir 05_OG.Inference_Phylogenomic/05_tree/00_species
conda activate tree
iqtree -m TESTNEW -b 100 -s conc_species_tree --prefix 05_tree/00_species/species_tree -nt 9
```

## Alignment and trimming on 02_disco_OG folder 
```
mkdir 03_aligned/01_prova_alignment_02_disco_OG
mkdir 04_trimmed/01_prova_trimmed_02_disco_OG
```
Alignment
```
cd  05_OG.Inference_Phylogenomic/02_disco_OG
conda activate sequence
for file in *.faa; do mafft --auto --anysymbol "$file" > ../03_aligned/01_prova_alignment_02_disco_OG/${file/.faa/_aligned.faa}; done
```

Trimming
```
cd ../03_aligned/01_prova_alignment_02_disco_OG
conda deactivate #you need to deactivate all the open environments to run 'bmge'
for file in *.faa; do bmge -i "$file" -t AA -m BLOSUM62 -e 0.5 -g 0.4 -of ../../04_trimmed/${file/_aligned.faa/_trimmed.faa}; done
```

