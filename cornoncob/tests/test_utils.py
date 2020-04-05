import csv
import sys
import random
import subprocess

sys.path.append(".")

import numpy as np

from Bio import SeqIO

from io_utils import parse_gff, if_not_exists_make

TEST_PEPS = './test_peps.csv'
# path should remain same for testing not user given


def read_test_peps():
    '''
    Reads in peptides and their consensus nucleotide translation
    into tuples and returns those tuples in a list. This is the ground truth
    list.
    '''
    with open(TEST_PEPS) as tp:
        reader = csv.reader(tp)
        next(reader)
        return reader.readrows()


def run_baseline_prokka(phenotype, run_dir):
    '''
    Runs prokka on all the genomes present in the given phenotype and then
    copies or stores in memory the results of each prokka prediction gff file
    in order to establish where the non-coding regions of the genomes will be
    predicted to be. This way test peptide coding nuceltide sequences can be
    safley inserted into each genome.
    '''
    test_files = if_not_exists_make(run_dir, 'test_files')
    for genome in phenotype:  # run prokka on all genomes in the phenotype
        genome.make_gene_predictions()
    
    # copy genome files into the test files folder to conserve them
    # after the run is complete maybe dont need to do this.
    
    # run prokka on each insert into the genomes with those files and
    # then just let them get overwritten. probably want to store some
    # info on what was inserted where so can check that out later on
    


def insert_test_peps(positive_pheno, test_peps, add_stop_codon=False):
    '''
    Given a phenotype which is considered the positive for a given character
    inserts the nucleotide sequences of the peptides in TEST_PEPS file into the
    non-coding regions of the genome as defined by the gene predicitions made by
    prokka and written to gff file.
    
    :param positive_pheno: Phenotype. Phenotype object to use as test case
    :param test_peps: List. List of lists of test_peps.csv file
    :param add_stop_codon: Boolean. If True inserts stop codon at start of test peptides
    '''

    num_test_peps = len(test_peps)
    for genome in positive_pheno:
        genome_records_dict = SeqIO.to_dict(
            SeqIO.parse(genome.genome_file, 'fasta'))
        gene_records = parse_gff(genome.gff_path, 0, 3, 4)
        random_indicies = np.random.choice(
            len(gene_records), num_test_peps, replace=False)

        random_gene_headers = [gene_records[i] for i in random_indicies]
        # pick random records predicted to be coding. Index 0 is header name
        # in the genome fasta 1 = start of coding, 2 = end coding

    for i, gene_header, start, end in enumerate(random_gene_headers):
        peptide = test_peps[i]
        # start and end refer to start end of coding region in the record
        # buffer is distance between test pep and coding region
        # where = random.randint(0, 1) adding start later
        buffer = random.randint(25, 30)  # probably want to make relative
        in_seq = genome_records_dict[gene_header].seq
        genome_records_dict[gene_header].seq = in_seq[:end] + \
            in_seq[end:end+buffer] + peptide[2] + in_seq[end+buffer:]
        # insert the test pep nucl sequence

    # overwrite the existing genome file with the modified and unmodified
    # sequence records
    SeqIO.write(genome_records_dict, genome.genome_file, 'fasta')
        
    
