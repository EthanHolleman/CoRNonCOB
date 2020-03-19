from genome import Genome

from io_utils import if_not_exists_make


ECOLI = '/homes/aahmad9/GCF_000005845.2_ASM584v2_genomic.fna'
OUT_DIR = '/homes/aahmad9'


new_genome = Genome(ECOLI, 'A')
print(new_genome)

new_genome.make_gene_predictions(OUT_DIR, 2, 'prokka','prokka_results')

