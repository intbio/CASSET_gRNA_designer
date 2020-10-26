# CASSET_gRNA_designer
A library for gRNA design (for both single Cas and split-systems).  
Split-system means 2 Cas (from the list: SpCas9, SaCas9, CjCas9, StCas9) are located in orientation (PAM-out, PAM-in, PAM-direct-forward, PAM-direct-reverse) with the condition: the angle between protein ends (N- or C-) less than 20 degrees.  
Also, this library allows to design gRNA to recognize targets in a set of pathogens.
Parameters that can be changed:
- permissible angle between 2 ends of Cas proteins for our set of the distance between 2 Cas
- manually set distance between 2 Cas
- number of mismatches in the first 12 nt and in the whole targets (for the task to detect several pathogens) for both left and right.

### Pipeline:
![alt text](https://github.com/intbio/CASSET_gRNA_designer/blob/main/pipeline.png)

### Examples
In [example notebook](example.ipynb) you can find the following examples:
- search targets in N and E genes in Sars-Cov-2 genome for all variants of split-systems and for CasX
- search targets for NS5A gene and for all HCV genome for all variants of split-systems and for CasX
- search targets to detect 4 species of Borrelia simultaneously for all variants of split-systems

### Under development
- On-target and off-targets scores for designed gRNAs, ranging of gRNAs
- BLAST check of organisms in which designed gRNAs could work

#### Attribution
1. Code - A. Gribkova
2. Distance between Cas proteins in split-systems from molecular modeling - R. Novikov, G. Armeev
3. Logo - G. Armeev
___________________


![alt text](https://github.com/intbio/CASSET_gRNA_designer/blob/main/logo.png)
