# Implementing, adapting and applying the transPACT software on the VlaamseSuperComputer (VSC)
<i> <b>Master's thesis - Laboratory for Biomolecular Discovery and Engineering</b>, MSc in Bioinformatics, VIB - KU Leuven, 2022, Antoine Ruzette, Joleen Masschelein</i>

All credits regarding the transPACT platform goes to the authors of the following repository: https://github.com/chevrm/transPACT.
<br><br>

**Goal of the project:** 
<hr style="border:2px solid blue">
<ul>
  
  - Use the transPACT freeware as an evolutionary comparison between trans-AT query KS sequences and KS sequences of reference
  - Functional annotation of the query KS sequences 
  - Identification of the phylogenic relations of the query KS sequences (through the generation of a dendrogram phylogeny)
  - Genome mining approach: identification of the phylogenetic relationships of the antiSMASH database i.e. about 3 000 000 domain sequences distributed across 1937 trans-AT polyketide bacterial cluster 
  <br><br>
  
</ul>

**Implementation:**
<hr style="border:2px solid blue">

<ul>
  A comprehensive tutorial is available <a href="https://antoineruzette.github.io/BDE-transPACT-implementation/docs/TransPACT%20tutorial.html">here</a>.
  <br><br>  
</ul>

**Dissertation:**
<hr style="border:2px solid blue">

<ul> 
  My master's thesis dissertation can be accessed <a href="https://antoineruzette.github.io/BDE-transPACT-implementation/docs/paper/Ruzette,%20Masschelein,%202022.pdf">here</a>.
  <br><br>
</ul>

<b> Results: </b>
<hr style="border:2px solid blue">

<ul> <b>Query Targeted Approach:</b> Phylogenetic relationships of GO485 (Massilia <i> flava </i>), AWB69 (Caballeronia <i> udeis </i>) and PBC (Pseudomonas <i> baetica </i>) with 647 KS sequences of reference <br><br>

Update on April 10, 2022: addition of DK441 (Aquimarina sp. AU58 NZ_OMPD010000014), XB90 (Burkholderia <i> pseudomallei </i> MSHR4299 scaffold2), DO63 (Burkholderia <i> pseudomallei </i> strain BDE scaffold2 NZ_KN150938), BLU84 (Lysobacter <i> enzymogenes </i>), ASC93 (Massilia sp. Root335 LMCY010000011 DSM 102448), CR152 (Massilia <i> violaceinigra </i> B2 NZ_CP024608 DSM 19531) and K350 (Sporocytophaga <i> myxococcoides </i>). 
<br>Update on April 22, 2022: addition of Leinamycin, Guangnanmycin, Weishanmycin and Largimycin. 
<br>Update on June 21, 2022: addition of K350. 

- April 10, 2022: rectangular and circular phylogenies: https://itol.embl.de/export/94224211181134871648805232 (available in .svg and .pdf format) 
- April 22, 2022: rectangular and circular phylogenies: https://itol.embl.de/export/62205120139364271650624817
- June 21, 2022: rectangular and circular phylogenies: https://itol.embl.de/export/8520111949424021655827199 <br><br><br>

<b>Genome Mining Approach:</b> Phylogenetic analysis of publicly available PKSs from the antiSMASH ClusterBlast database, including GO485, AWB69 and related PKSs 
<br><br>
In order to reduce computing time and to concentrate efforts on C. <i> udeis </i> and M. <i> flava </i> related PKSs, only the 10 to 16 KS domains long PKSs have been selected from the antiSMASH ClusterBlast database. As a result, 511 trans-AT PKSs have been identified. The dendrogram visualizations 519 PKSs, including the 8 under study, can be accessed through an ITOL sharing link. 

- July 19, 2022: rectangular and circular phylogenies: https://itol.embl.de/export/85201185207360831658172247 
</ul>
