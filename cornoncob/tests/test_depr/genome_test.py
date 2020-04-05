'''
import os
import sys
import inspect
sys.path.insert(1, os.path.realpath(os.path.pardir))

from genome import Genome


LC40_GENOME = '/home/ethan/Documents/ecoli_genome/putonti_seqs/killers/Lc20.fasta'
LCX_GENOME = '/home/ethan/Documents/ecoli_genome/putonti_seqs/nice/Lc9873.fasta'
PROKKA = '/home/ethan/prokka/bin/prokka'


forty_g = '/home/ethan/Documents/github/CoRNonCOB/tests/Lc40.fasta/PROKKA_03212020.gff'
x_g = '/home/ethan/Documents/github/CoRNonCOB/tests/Lc9873.fasta/PROKKA_03212020.gff'

genomes = [
    Genome(LC40_GENOME, 'nice'),
    Genome(LCX_GENOME, 'nice')
]

genomes[0].gene_prediction_file = forty_g
genomes[1].gene_prediction_file = x_g

for g in genomes:
    g.get_non_coding_regions()
    g.translate_non_coding_seqs()
    g.write_peptides_to_fasta_file('.')

    cmd = ['cd-hit-est', '-i', input_file, '-o', output,
           '-c', str(identity),  '-T', '6', '-d', '0']
# cd hit guide command for keeping headers in tact
'''