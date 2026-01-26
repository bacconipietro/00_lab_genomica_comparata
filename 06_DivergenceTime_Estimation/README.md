# Divergence Time Estimation

Pipiline keeps going with ultrametric-tree task. To complete next passages we need to build a time tree from our phylogenomic inference. Traditionally, a Bayesan analysis using **BEAST** program is recommended, but in this course divergence dating is a secondary goal to reach another result. For this reason we use IQTREE with Least Square Dating (LSD2) algorithm wich uses a linear regression approach to fit a molecular clock to the branch lengths of a phylogenetic tree. It reconciles the observed substitution rates with user-provided date constraints, such as fossil calibrations or tip date,s to estimate the divergence times of all nodes.


In our case we will use Ancestral Dating to calibrate our tree using divergence times retrieved from public databases. We obtain tips dates data uploading our dataset on [TimeTree.org](https://timetree.org/). Make sure to upload only species name without undescores. 

```bash
cut -f2 dataset.tsv | sed ‘s/_/ /’
```

Next it's made a `calibration.txt` file where are reported tips dates and taxon involved, using only relevant nodes with this structure:

```
taxon1,taxon2 -50
taxon3,taxon4,taxon5 -100
```

Then IQ-tree starts:

```bash
conda activate time
#[time]
ln -s /home/STUDENTI/pietro.bacconi/Lab_genomica_comparata/00_practice/05_OG.Inference_Phylogenomic/04_trimmed/00_single_complete/conc_species_tree
iqtree -s conc_species_tree --date calibration.txt --date-tip 0 -o Toxrut -m Q.INSECT+F+I+R3 -nt 13 --prefix time_tree --date-options "-u 1"
```
