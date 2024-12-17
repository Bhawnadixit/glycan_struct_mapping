# glycan_struct_mapping
This is an open-source repository for converting specific glycan IUPAC strings into structure mapping for computing phi/psi (φ/Ψ) dihedrals. Thw data in this repo is taken from Reference:

Dixit, B., Vranken, W., & Ghysels, A. (2024). Conformational dynamics of α‐1 acid glycoprotein (AGP) in cancer: A comparative study of glycosylated and unglycosylated AGP. Proteins: Structure, Function, and Bioinformatics, 92(2), 246-264.

## Introduction

This document explains how to convert an **IUPAC glycan sequence** into a **3D structure mapping** of the glycan chain. The example glycan chain is as follows:

> Gal(b1-3)GlcNAc(b1-2)Man(a1-3)[Gal(b1-3)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-ASN15

This sequence represents a **branched glycan**, where the core sugars are linked together with specific **glycosidic linkages**, and branching occurs at certain points in the chain.

For standard nomenclature of glycans, please refer to Symbol Nomenclature for Glycans (SNFG)(https://www.ncbi.nlm.nih.gov/glycans/snfg.html).


#### 𝛃Gal(1→3)-βGlcNAc(1→2)-αMan(1→6)  (N-glycan branch1)
#### -------------------------------|           (N-glycan core)
#### -----------------------------𝛽Man(1→4)-𝛽GlcNAc(1→4)-𝛽GlcNAc(1→)-ASN15  
#### -------------------------------|  
#### βGal(1→3)-βGlcNAc(1→2)-αMan(1→3)  (N-glycan branch1)
The 2D representation of a glycan is shown as above, here α/β represent the isomers of the monosaccharide, 1→3 represents glycosidic linkage between two monosaccharides, ASN15 is the asparagine 15 amino acid residue of the protein α-1 acid glycoprotein (AGP).

Here, α refers to the configuration of a cyclic sugar in which the oxygen on the anomeric carbon is positioned on the opposite side of the ring compared to the substituent on the adjacent carbon. This is in contrast to β, where the oxygen and the adjacent substituent are on the same side of the ring.

In N-glycans, monosaccharides can be connected through various glycosidic bonds, such as 1→1, 1→2, 1→3, 1→4, 2→3, 1→6, or 2→6, where the numbers refer to the carbon atoms of the two involved monosaccharides. The relative orientations of the monosaccharides are defined by torsion angles: ϕ and ψ for 1→1, 1→2, 1→3, 1→4, and 2→3 linkages, and ϕ, ψ, and ω for 1→6 and 2→6 linkages. For instance, a disaccharide formed by α-Neu5A and β-Gal linked through a 2→6 glycosidic bond. For example, for α-Neu5A-(2→6)-β-Gal, the three torsion angles are—ϕ (O6-C2-O6-C6), ψ (C2-O6-C6-C5), and ω (O6-C6-C5-O5).

## Problem
N-glycans can exhibit varying degrees of branching, such as bi-antennary, tri-antennary, or tetra-antennary structures. This branching leads to different indexing patterns; for example, starting from the core, the first three monosaccharides are numbered 1, 2, and 3, with branch 1 continuing to 4, 5, and 6, while branch 2 starts from 7, 8, and 9. Depending on these branching indices, the selection of φ (phi) and ψ (psi) angles in the glycan chain can become complex. Accurate information about the glycosidic linkages between the saccharides is crucial to ensure proper structure and to avoid manual errors in selecting these linkages in post-Molecular dynamics analysis.

Therefore, this repository converts the IUPAC name into triplets, to identify branches, glycosidic angles, and appropriate indices, and then converts each monosaccharide-monosachharide glycosidic linkage, residue names (according to GROMACS/PyMol based naming conventions, indices into a dataframe.

## For the full demo please see the Jupyter notebook, test_glycans.ipynb

## Python files in ./iupac_to_mapping/

- **string_process.py**: python script to convert IUPAC string into a list of triplets
- **iupac_converter.py**: python script to convert triplets into core and branch structure.
- **glycan_chain_indices.py**: python script to identify and list all glycans in IUPAC string and map them in the .pdb file
- **compute_torsions.py**: python script to compute  ϕ, ψ, and ω angles either on a single .pdb, or a trajectory file.

## Requirements
Using 'requirements.txt' or 'environment.yml'
> pip install requirements.txt

> conda env create -f environment.yml

> conda activate glycan_env

**NOTE**: This is not tested on other glycans, therefore, to run the code can be adapted. 


