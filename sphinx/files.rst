Output File Structure
=====================

This is what I am thinking the output file structure should look like ish

.. code-block:: text

    .
    └── run_name
        ├── conserved_motifs.fasta
        ├── phenotype_a
        │   ├── genome_a
        │   │   ├── non_coding_seqs.fasta
        │   │   └── prokka_results
        │   ├── genome_b
        │   └── genome_c
        ├── phenotype_b
        │   ├── genome_d
        │   ├── genome_e
        │   └── genome_f
        └── unique_motifs.fasta

User supplies the parent directory and the run name and the program creates
the rest of the directory structure which will approximately reflect the way
objects are set up in the pipeline. 

conserved_motifs.fasta and unique_motifs.fasta just represent the final
output that the user would be interested in. Genome_a is being used as the
example genome; all genomes would have similar structure with a sub directory
holding their prokka results and a file with the non coding regions.
