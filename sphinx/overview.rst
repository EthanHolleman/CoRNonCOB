Project Overview
================


Abstract
---------
CoRNonCOB (Compare Really Non-Conserved Oligopeptides in Bacteria) is a tool
that identifies peptides in the non-coding regions of bacterial genomes
that may contribute to a phenotype of interest. CoRNonCOB is designed to work
at the very specific scale of individual phenotypes of specific strains of
bacteria and because of this RNA-seq data that could illuminate the expression
profiles of the bacteria of interest may not be available. Therefore,
CoRNonCOB attempts to leverage genomic data by identifying peptides in
non-coding regions that are conserved in the members of a positive phenotype
and absent in a control.

With this tool, we hope to identify the peptides that are unique
to *Lactobacillus crispatus* that produce a microbial
agent to kill E. coli bacteria. Through identifying these unique peptides,
we hope to provide another parameter to characterize phenotypes. 
Our tool takes x amount of genomes in FASTA format  from each
phenotype, and generates the peptides unique to each 
respective phenotype. Our tool is available at the
`GitHib page here <https://github.com/EthanHolleman/CoRNonCOB>`_.

For More Information
--------------------

For more information on the details of the project, pipeline structure,
results and citations please visit our application note which can be downloaded
from the link below.

CoRNonCOB Application Note :download:`pdf <../graphs_and_misc/app_note/App_Note.pdf>`