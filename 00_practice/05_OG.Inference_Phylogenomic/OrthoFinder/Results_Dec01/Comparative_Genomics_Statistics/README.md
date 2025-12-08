## Move Orthofinder results folders
```
mv Orthofinder ../../05_OG.Inference_Phylogenomic
mkdir 05_OG.Inference_Phylogenomic/01_DISCO
```

# Paralogy filtering 
## Run Disco
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
