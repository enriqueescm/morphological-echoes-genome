# morphological-echoes-genome
From Typology to Genomics: Revisiting Human Morphological Diversity Through Population Genetics
## Overview

Most bioinformatics projects treat human populations as statistical clusters.
This one asks what the genome reveals about the morphological diversity
that physical anthropologists have documented for over a century.

This project integrates **population genomics** with **biological anthropology**
to investigate whether craniofacial variation across human populations
leaves detectable signatures in the genome — through population structure,
natural selection, and candidate gene analysis.

## Biological Question

Can we identify genomic signals that parallel the craniofacial morphological
variation observed across human populations?

## Approach

| Module | Analysis | Tools |
|--------|----------|-------|
| 01 | Population structure (PCA, FST) | PLINK2, R |
| 02 | Selection scans (iHS, XP-EHH) | Python, selscan |
| 03 | Candidate gene analysis | R, Bioconductor |
| 04 | Final report | Quarto |

## Data Sources

- [1000 Genomes Project Phase 3](https://www.internationalgenome.org) — 2,504 individuals, 26 populations
- [GWAS Catalog](https://www.ebi.ac.uk/gwas) — craniofacial trait associations
- [UCSC Genome Browser](https://genome.ucsc.edu) — functional annotation

## Key Candidate Genes

`EDAR` `PAX3` `TP63` `BMP4` `ALX1` `RUNX2` `MSX1` `COL1A1`

These genes have documented roles in craniofacial development and show
population-differentiated allele frequencies consistent with local adaptation.

## Author

**Enrique Estevez**
Biological Anthropologist & Biomedical Researcher  
MSc Bioinformatics & Biostatistics (in progress)  
Interests: population genomics, morphological evolution, RWE
[Email](mailto:enrique.escm@gmail.com)

## Status

🟡 In progress — Module 01 (Population Structure) active
