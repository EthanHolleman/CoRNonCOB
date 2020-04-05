from json import JSONEncoder

GFF = '/home/ethan/Documents/github/CoRNonCOB/corncob/killers/Lc20.fasta/prokka_results/PROKKA_03222020.gff'

GENOMES = '/home/ethan/Documents/ecoli_genome/putonti_seqs/nice'
RUN_DIR = '/home/ethan/Documents/phenotype_test'
PROKA = '/home/ethan/prokka/bin/./prokka'

from phenotype import Phenotype

p = Phenotype(GENOMES, RUN_DIR, phenotype='n')
p.pull_peptides(prokka_exec=PROKA)
p.get_conserved_sequences()

