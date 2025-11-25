# Working on Dataset 

## Download dataset 
```
cd 00_pratice/00_data
conda activate sequence
bash ../99_scripts/download_dataset.sh dataset.tsv
```

## AGAT process
#### Keeping longest isoform
```
conda activate GAAS
cd 01_gff
for gff in *.gff; do
agat_sp_keep_longest_isoform.pl -gff "$gff" -o ${gff/.gff/_longest.gff}
done 
```
#### Extracting and translating sequences 
```
cd 00_data
mkdir 02_raw_proteoms
for gff in 01_gff/*_longest.gff; do
agat_sp_extract_sequences.pl -g "$gff" -f ../00_genome/${gff/_longest.gff/.fna} -t cds -p --cfs --output 02_raw_proteoms/${gff/_longest.gff/faa}
done
```

## Remove pseudogenes
```
cd 02_raw_proteoms
bash ../../99_scripts/pseudogenes_find_eliminate.sh
for pseudo in *.txt, do wc -l "$pseudo" done
```

## Modify headers
```
# [inside 02_raw_proteoms]
for prot in *.faa; do
ID=$(basename -s .faa "$prot")
sed -i.old -E "s/>(rna-XM_[0-9]+\.[0-9]) (gene=gene-.[^ ]+) name=(.[^ ]+) .+$/>${ID}\|\3/" "$prot"
done

rm raw_proteoms
rm *.old
```
