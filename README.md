# Alignement_construction

>**Authors**: Sophie Gallina 1, Zoé Postel 1 <br />

>**Affiliation**: CNRS, Univ. Lille, UMR 8198 – Evo-Eco-Paleo, F-59000 Lille, France <br />

&nbsp;

Biopyhton scripts developed with Sophie Gallina for the construction of multifasta files from a vcf file and a reference sequence. Briefly, the script takes as input a reference sequence in fasta format + a variant file in vcf format. For each individual and each gene identified in the vcf file, the script will take the reference sequence and modify it by writing the correct nucleotide for each position identified as variable in the vcf file. This sript was initially elaborrated to reconstruct mitochondrial gene sequences (supposedly haploid) (Postel et al, in prep). 

Three scripts are given, which mainly does the same purpose but with slight variation to deal with the heterozygote position (i.e. more than one allele at a position): <br />
-- Alignement_construction.py : version to mask the heterozygote position with a "N" <br />
-- Alignement_construction_allele1.py : version to write, at the heterozygote positions, the allele 1 <br />
-- Alignement_construction_allele2.py : version to write, at the heterozygote positions, the allele 2 <br />

&nbsp;

_Command line to run the script:_

./script_name.py --name xx  --AD_thr xxx

> --name: give the name of the analyse (i.e. name of the folder containing the reference fasta file and the vcf + the name of the fasta file).  <br />
> --AD_thr: minimum number of read per allele to modify the fasta sequence with the main variants. 


