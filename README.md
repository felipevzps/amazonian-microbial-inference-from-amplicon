# sample-inference-from-amplicon

Sample inference from amplicon data with single-nucleotide resolution using DADA2.

## Dataset
Six soil samples collected from the Amazonian forest and pasture were sequenced using the Illumina MiSeq platform in PE 250bp mode. The sequencing employed the primers 515F (Parada et al. 2016) and 806R (Apprill et al. 2015), with lengths of 19 and 20 bp, respectively.
[Earth Microbiome Protocols](https://earthmicrobiome.org/protocols-and-standards/16s/)

## Amplicon Sequence Variants (ASVs)
| Samples | ASVs  |
|---------|-------|
| 6       | 2179  |
>The distribution of sequence lengths shows the presence of 2,179 ASVs.
>
>The majority having a length of 253bp (1,954 ASVs).

## ASVs after chimera removal
| Samples | ASVs  |
|---------|-------|
| 6       | 2125  |
>Chimera / non-chimera ratio: 98.1663%
>
>Chimera sequences represent only ~2% (54) of the merged sequences.

## Sequence control by each step
| Sample | Input | Filtered | DenoisedF | DenoisedR | Merged | Nonchim |
|--------|-------|----------|-----------|-----------|--------|---------|
| A1     | 108005| 55786    | 54220     | 53897     | 45782  | 44818   |
| A2     | 90233 | 47722    | 46110     | 46064     | 38333  | 37368   |
| A3     | 117013| 60818    | 59185     | 59044     | 50321  | 49122   |
| A7     | 133600| 68871    | 65140     | 65104     | 48824  | 48070   |
| A8     | 60867 | 31178    | 28789     | 28872     | 21239  | 21029   |
| A9     | 93554 | 46895    | 43851     | 43844     | 33108  | 32843   |

>This show a satisfactory progression, with the primary loss of sequences occurring during raw read filtering for quality and size.
>
>However, a substantial portion of the clean reads is retained in subsequent pipeline steps.

## Accessing ASV Taxonomy

The `assignTaxonomy` function classifies sequences using a reference set of sequences with known taxonomy. In this project, the Silva database `silva_nr99_v138.1_train_set.fa` was used as reference.
>Download can be done [here](https://zenodo.org/record/4587955)

```R
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa", multithread=FALSE)
write.table(taxa, "results/taxa.txt", sep="\t", quote=F)

# Extract sequences for easy visualization
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

#       Kingdom    Phylum           Class                 Order               Family               Genus
# [1,] "Bacteria" "Proteobacteria" "Alphaproteobacteria" "Rhizobiales"       "Xanthobacteraceae"  NA   
# [2,] "Bacteria" "Proteobacteria" "Alphaproteobacteria" "Rhizobiales"       "Xanthobacteraceae"  NA   
# [3,] "Archaea"  "Crenarchaeota"  "Nitrososphaeria"     "Nitrososphaerales" "Nitrososphaeraceae" NA   
# [4,] "Bacteria" "Proteobacteria" "Alphaproteobacteria" "Rhizobiales"       "Xanthobacteraceae"  NA   
# [5,] "Archaea"  "Crenarchaeota"  "Nitrososphaeria"     "Nitrososphaerales" "Nitrososphaeraceae" NA   
# [6,] "Archaea"  "Crenarchaeota"  "Nitrososphaeria"     "Nitrososphaerales" "Nitrososphaeraceae" NA
```
>Proteobacteria and Crenarchaeota are the most abundant phylum.

## Rarefaction curve after sample filtering using truncQ=2

<img src="https://github.com/felipevzps/sample-inference-from-amplicon/blob/main/results/rarefaction_Q2.PNG" width="800">

>With default settings (truncQ=2), samples A7 and A9 exhibit the highest species diversity, with over 1000 and 700 species, respectively.

## Rarefaction curve after sample filtering using truncQ=20

<img src="https://github.com/felipevzps/sample-inference-from-amplicon/blob/main/results/rarefaction_Q20.PNG" width="800">

>Increasing the quality filtering parameter (truncQ=20) becomes more restrictive, discarding over 96% of initial reads and resulting in a significant reduction in detected species.

## Alpha diversity analyzed by Microbiome Analyst

<img src="https://github.com/felipevzps/sample-inference-from-amplicon/blob/main/MicrobiomeAnalyst/3_alpha_diver_8.png" width="516"/> <img src="https://github.com/felipevzps/sample-inference-from-amplicon/blob/main/MicrobiomeAnalyst/4_alpha_diverbox_8.png" width="484"/>

>The alpha diversity graph at the Phylum level indicates two distinct groups (A and B) based on similar species composition within each group.

## Taxonomy of alpha diversity

<img src="https://github.com/felipevzps/sample-inference-from-amplicon/blob/main/MicrobiomeAnalyst/2_taxa_alpha_5.png" width="800">

>Examining alpha diversity taxonomy reveals that groups A and B, while containing similar species, exhibit significantly different proportions of these species.

## Beta diversity analyzed by Microbiome Analyst

<img src="https://github.com/felipevzps/sample-inference-from-amplicon/blob/main/MicrobiomeAnalyst/1_PCoA_3D_beta_diversity.png" width="800">

>Beta diversity analysis confirms the differences between groups A and B, indicating distinct diversity at the Phylum level.

## Treatment Bioindicators = LEfSe (LDA score)

<img src="https://github.com/felipevzps/sample-inference-from-amplicon/blob/main/MicrobiomeAnalyst/5_bar_graph_6_LEfSe_phylum_relevant_treatments.png" width="800">

>Due to differences in microbiome composition, groups A and B also show different functional diversity, with group A associated with nitrogen fixation and nutrient cycling, while group B is linked to methane cycling and cellulose/sulfate reduction.

## Observations

- All samples consistently reached the plateau in the rarefaction curve (default - truncQ=2) , indicating robust sequencing coverage.
- Stringent high-quality sequencing filtering criteria (e.g., Q>20) resulted in a reduced number of sequences and fewer unique ASV identifications. Rigorous quality filtering revealed around 200 ASVs in the rarefaction curve, emphasizing the trade-off between accuracy and the ability to detect rare taxa.
- Taxonomic annotation of ASVs allowed the detection of alpha and beta diversities, revealing differences between treatments A and B in terms of microbial community taxa and their distinct functions in Amazonian forest and pasture.
