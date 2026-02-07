# Dataset creation

To conduct a comparative genomic analysis we need genomes species to compare. We start creating a dataset with our model organism and 4/5 other species, chosing them based on our biological study purpose.
Here it's presented a mosquitos genomes [dataset.tsv](https://github.com/bacconipietro/00_lab_genomica_comparata/blob/main/00_practice/00_data/dataset.tsv) searched on NCBI, the model we chose to study is [**_Anopheles stephensi_**](https://mesamalaria.org/resource-hub/pmi-anopheles-stephensi-resource-page/). 

Mosquitoes are predominantly tropical and subtropical vectors, yet several lineages have successfully colonized temperate regions. This ecological transition requires for sure significant physiological and genomic modifications to survive overwintering, freeze tolerance, and altered metabolic rates. To identify the genomic signatures we selected two species known to inhabit higher latitudes or cooler altitudes, and three one restricted to equatorial or warm regions. 

| AN | Species | Code | Regions |
| :---: | :---:	| :---: | :---: |  
| GCA_000441895.2 | *Anopheles_sinensis* | Anosin | temperate |
| GCF_029784135.1 | *Toxorhynchites_rutilus_septentrionalis* | Toxrut | temperate |
| GCF_943734745.1 | *Anopheles_darlingi* | Anodar | tropical |
| GCF_943734735.2 | *Anopheles_gambiae* | Anogam | troipcal |
| GCF_943734845.2 | *Anopheles_funestus* | Anofun | tropical | 





-----

## Download
We download genomes data from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=7164&assembly_level=1:3), collecting accession numbers we set some filters:

+ Reference genome
+ Annotated genome
+ Scaffold as minimum level for the assembly
+ Released after 2010

Here is reported the script: 
```bash
#!/bin/bash

#This script downloads genomes and relative gff from NCBI datasets creating ready to use folders

AN2name=$1

mkdir 00_genome
mkdir 01_gff

while IFS=$'\t' read -r AN sname ID; do
	echo $AN
	#download specifying the name of the folder
	datasets download genome accession "$AN" --filename "$ID".zip --include genome,gff3
	#unzip specifying the name of the folder
	unzip "$ID".zip -d "$ID"
	#rename the two file of interest
	mv "$ID"/ncbi_dataset/data/"$AN"/*.fna 00_genome/"$ID".fna
	mv "$ID"/ncbi_dataset/data/"$AN"/*.gff 01_gff/"$ID".gff
	#delete the folder
	rm -r "$ID"/
done < "$AN2name"
```
**Run download**
```bash 
nano dataset.tsv
conda activate sequence
bash download_dataset.sh dataset.tsv
```
-----

## Keeping longest isoform  

Following genome annotation, it is common to find multiple isoforms for a single gene, especially when transcriptomic data is used to improve precision. However, retaining all isoforms can bias downstream comparative analyses. To resolve this, we use the AGAT suite (Another GFF Analysis Toolkit) to filter the annotation. We run the Perl script `agat_sp_keep_longest_isoform.pl` to discard shorter variants, retaining only the single longest representative for each gene. Once filtered, we use `agat_sp_extract_sequences.pl` to parse the GFF file, extracting the corresponding nucleotide sequences and translating them into protein sequences for functional analysis.

```bash
#[GAAS]
cd 01_gff
for gff in *.gff; do
agat_sp_keep_longest_isoform.pl -gff "$gff" -o ${gff/.gff/_longest.gff}
done 
```

```bash
#[GAAS]
cd 00_data
mkdir 02_raw_proteoms
for gff in 01_gff/*_longest.gff; do
agat_sp_extract_sequences.pl -g "$gff" -f ../00_genome/${gff/_longest.gff/.fna} -t cds -p --cfs --output 02_raw_proteoms/${gff/_longest.gff/faa}
done
```
-----

## Remove pseudogenes

To ensure our dataset contains only functional proteins, we must detect and remove sequences labeled as coding regions that contain internal STOP codons, **pseudogenes**.
Here it's used a full script to automate this cleaning process by discarding any sequence containing a stop codon. 

Here is the script:
```bash
#!/bin/bash/

# list each plausible pseudogene present. 

mkdir raw_proteomes
mv *.faa raw_proteomes/

#make proteome one line
cd raw_proteomes
for proteome in *.faa; do
	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' "$proteome" > ../${proteome/.faa}".faa"
done
cd ..

#Extract all pseudogene names
mkdir 00_pseudogene_name

for proteome in *.faa; do
	species=$(basename -s .faa "$proteome")
	grep -B1 '*' "$proteome" | grep ">" >> 00_pseudogene_name/"$species"_pseudogenes_name.txt
done

#removes sequences identified as pseudogenes

for pseudo_file in 00_pseudogene_name/*_pseudogenes_name.txt; do
	species=$(basename -s _pseudogenes_name.txt "$pseudo_file")
	while IFS=$'\t' read -r header; do
		sed -E -i "/${header}/{N;d;}" "$species".faa # N option loads the next line found after the pattern and put it into pattern space too; d delete the pattern space
	done < "$pseudo_file" 
done

mv 00_pseudogene_name ../00_genome
```

**Then run:**
```bash
cd 02_raw_proteoms
bash ../../99_scripts/pseudogenes_find_eliminate.sh
for pseudo in *.txt, do wc -l "$pseudo" done
```
-----

## Modify headers

We will modify the headers to retain only the species identifier. This simplification is indispensable, ensuring that subsequent software can correctly track and identify each gene without ambiguity.

```bash
for prot in *.faa; do
ID=$(basename -s .faa "$prot")
sed -i.old -E "s/>(rna-XM_[0-9]+\.[0-9]) (gene=gene-.[^ ]+) name=(.[^ ]+) .+$/>${ID}\|\3/" "$prot"
done    #Unfortunatly the pattern for Anosin.faa headers is different from 'rna-XM', for this reason the sed didn't run correctly. I will correct the single tip in the final tree file.
```
