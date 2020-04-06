Understanding the Output File Structure
=====================

Tree Overview
-------------
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

What is Significant to the User?
--------------------------------
Files that the user will be most interested in will be found in the first
layer of their run directory. These will include fasta files of candidate
peptides of interest, csv type file to piece together where everything came
from and a log file to help out with that and for debugging.

All files that are stored at deeper level is stuff Corn on Cob is using to do the actual
analysis and would include output from Prokka, other software used and temp
files. We are planning on leaving these intact for the final version for easier
debugging and in case the user wants to review and of the intermediate files.

