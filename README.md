[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1231518.svg)](https://doi.org/10.5281/zenodo.1231518)

# Pom152

This repository contains the modeling files and the analysis related to the
article ["Molecular Architecture of the Major Membrane Ring Component of the Nuclear Pore Complex"](https://www.ncbi.nlm.nih.gov/pubmed/28162953)
by Upla et al. in Structure 2017. The scripts work with the [IMP](http://salilab.org/imp) (version 2.6).

## List of files and directories:

- `data`		            contains all relevant data 
   
   `*.pdb` : Representation PDB for individual domains

   `*.csv` : Cross-links data for further validation

   `FASTA_*.txt` : Sequence files

  `data/Cadherin` : sequence and structure alignments between Cadherin and Pom152
  
  `data/em2d_T1_truncation_1-1135` : EM 2D class averages of the Pom152 1-1135 (truncated)
  
  `data/em2d_T2_truncation_1-936` : EM 2D class averages of the Pom152 1-936 (truncated)
  
  `data/em2d_WT` : EM 2D class averages of the full length Pom152
  
  `data/em3d_density` : EM 3D density map of the full length Pom152
  
  `data/HHPred` : Template search results using HHPred
  
  `data/MODELLER_all` : Comparative models generated using MODELLER
  
  `data/Pom152_comparative_models_by_SAXS` : Best-scoring comparative models filtered by the corresponding SAXS data


- `template`			                  contains modeling scripts

  modeling_pdb375_482.py  : modeling script


- `results`		                      contains resulting structures and output files

  `results/EM2D_selected` : Validation of the Pom152 structure using selected EM 2D class averages (using multiple resolutions ranging 35 - 80 Å)
  
  `results/Pom152_em3d_Final` : Best-scoring Pom152 structures (clusters 0 and 1) and localization probability density maps
  
  `results/Pom152_em3d_with_Rotational_Restraint` : Refined Pom152 structure using a dihedral restraint (holding a 90 degree rotation between neighboring Ig-fold domains)


- `prefilter`			                 Selection of 500 good-scoring models ranked by the combined total score.  Also score log files (*.log) are included

- `analysis`			                  contains Clustering scripts and results


## Information

_Author(s)_: Seung Joong Kim

_Date_: March 2017

_License_: [CC-BY-SA-4.0](https://creativecommons.org/licenses/by-sa/4.0/legalcode).
This work is freely available under the terms of the Creative Commons
Attribution-ShareAlike 4.0 International License.

_Last known good IMP version_: [![build info](https://integrativemodeling.org/systems/?sysstat=27&branch=master)](https://integrativemodeling.org/systems/) [![build info](https://integrativemodeling.org/systems/?sysstat=27&branch=develop)](https://integrativemodeling.org/systems/)

_Publications_:
 - Upla P, Kim SJ, Sampathkumar P, Dutta K, Cahill SM, Chemmama IE, Williams R, Bonanno JB, Rice WJ, Stokes DL, Cowburn D, Almo SC, Sali A, Rout MP, Fernandez-Martinez J. [Molecular Architecture of the Major Membrane Ring Component of the Nuclear Pore Complex](https://www.ncbi.nlm.nih.gov/pubmed/28162953), Structure, 2017, 10.1016/j.str.2017.01.006.
