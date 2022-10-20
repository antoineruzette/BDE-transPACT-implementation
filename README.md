# Implementing, adapting and applying the transPACT software on the VlaamseSuperComputer (VSC)
<i> <b>Master's thesis - Laboratory for Biomolecular Discovery and Engineering</b>, MSc in Bioinformatics, VIB - KU Leuven, 2022, Antoine Ruzette, Joleen Masschelein</i>

<b> Note that this is a private repository. Only persons who have been granted access should be reading this. </b> 

All credits regarding the transPACT platform goes to the authors of the following repository: https://github.com/chevrm/transPACT. 

**Goal of the project:** 
- Use the transPACT freeware as an evolutionary comparison between query KS sequences and KS sequences of reference
- Functional annotation of the query KS sequences 
- Identification of the phylogenic relations of the query KS sequences (through the generation of a dendrogram phylogeny)
- Genome mining approach: phylogenetic assessment of the 2000 PKS that are publicly available on the antiSMASH database


<b> Results: </b> 

**Query Targeted Approach:** Phylogenetic relationships of GO485 (Massilia flava), AWB69 (Caballeronia udeis) and PBC (Pseudomonas baetica) with 647 KS sequences of reference

Update on April 10, 2022: addition of DK441 (Aquimarina sp. AU58 NZ_OMPD010000014), XB90 (Burkholderia pseudomallei MSHR4299 scaffold2), DO63 (Burkholderia pseudomallei strain BDE scaffold2 NZ_KN150938), BLU84 (Lysobacter enzymogenes), ASC93 (Massilia sp. Root335 LMCY010000011 DSM 102448), CR152 (Massilia violaceinigra B2 NZ_CP024608 DSM 19531) and K350 (Sporocytophaga myxococcoides). 
Update on April 22, 2022: addition of Leinamycin, Guangnanmycin, Weishanmycin and Largimycin. 
Update on June 21, 2022: addition of K350. 

- April 10, 2022: rectangular and circular phylogenies: https://itol.embl.de/export/94224211181134871648805232 (available in .svg and .pdf format) 
- April 22, 2022: rectangular and circular phylogenies: https://itol.embl.de/export/62205120139364271650624817
- June 21, 2022: rectangular and circular phylogenies: https://itol.embl.de/export/8520111949424021655827199

**Genome Mining Approach:** Phylogenetic analysis of publicly available PKSs from the antiSMASH ClusterBlast database, including GO485, AWB69 and related PKSs 

In order to reduce computing time and to concentrate efforts on C. udeis and M. flava related PKSs, only the 10 to 16 KS domains long PKSs have been selected from the antiSMASH ClusterBlast database. As a result, 511 trans-AT PKSs have been identified. The dendrogram visualizations 519 PKSs, including the 8 under study, can be accessed through an ITOL sharing link. 

- July 19, 2022: rectangular and circular phylogenies: https://itol.embl.de/export/85201185207360831658172247 
