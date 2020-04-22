from cornoncob.io_utils import parse_gff, if_not_exists_make, convert_genome_to_header_dict
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import csv
import sys
import subprocess
import os

from cornoncob import TEST_PEPS
sys.path.append(".")


def read_test_peps(TEST_PEPS=TEST_PEPS):
    '''
    Reads in peptides and their consensus nucleotide translation
    into tuples and returns those tuples in a list. This is the ground truth
    list. Peptides with type P are positive and only inserted into the
    positive phenotype. Peptides with type C are control or negative and insert
    -ed into both genomes. We expect to see positive peptides in the final
    output but do not expect to see control / negative peptides.

    Returns a tuple of two lists, first is the positive peptides second is the
    control peptides.
    '''
    positive, control = [], []
    with open(TEST_PEPS) as tp:
        reader = csv.reader(tp)
        next(reader)  # skip header
        for row in reader:
            if row[-2] == 'C':
                positive.append(row)
                control.append(row)
            else:
                positive.append(row)
    return positive, control


def insert_test_peptides_into_all_phenotypes(postive_pheno, negative_pheno,
                                             seed=100, prokka_exec='prokka',
                                             stop=True):
    '''
    Given both the positive and negative phenotypes this function inserts
    the test peptide nucleotide sequences into the corresponding genomes.
    Positive peptides are only inserted into the positive phenotype genomes
    while negative peptides are inserted into all genomes.
    '''
    pos_test_peps, neg_test_peps = read_test_peps(TEST_PEPS)

    positive_test_data = insert_test_peptides_into_phenotype(
        postive_pheno, pos_test_peps, stop, seed, prokka_exec)
    negative_test_data = insert_test_peptides_into_phenotype(
        negative_pheno, neg_test_peps, stop, seed, prokka_exec)

    return positive_test_data, negative_test_data


def insert_test_peptides_into_phenotype(phenotype, test_peptide_list, stop=True,
                                        seed=100, prokka_exec='prokka'):
    '''
    Inserts the given test_peptides from the test_peptide_list into the
    given phenotype at random contigs. Returns a dictionary of contig headers
    and locations of peptide insertions.
    '''
    np.random.seed(seed)
    test_data = {}
    num_test_peps = len(test_peptide_list)

    for i, genome in enumerate(phenotype.genomes):
        genome_copy = f'{genome.genome_file}.copy'
        if '.copy' in genome.genome_file:
            continue  # copy of a genome do not add peptides
        else:
            genome.make_gene_predictions(path_to_exec=prokka_exec)
            genome_records_dict = SeqIO.to_dict(
                SeqIO.parse(genome.genome_file, 'fasta')
            )  # turn genome file into dictionary by header
            gene_records = parse_gff(genome.gene_prediction_file, 0, 3, 4)
            # pull out gene predictions
            random_indicies = np.random.choice(
                len(gene_records), num_test_peps, replace=False
            )  # pick random gene headers to insert into
            random_gene_headers = [gene_records[i] for i in random_indicies]
            # store the headers coresponding to random indicies

            for j, random_gene in enumerate(random_gene_headers):
                gene_header, start, end = random_gene[0], int(
                    random_gene[1]), int(random_gene[2])
                peptide = test_peptide_list[j]
                buffer = np.random.randint(25, 30)
                # num nucleotides between coding and inserted peptide start
                current_seq = genome_records_dict[gene_header].seq
                if stop:
                    peptide[2] = 'TAA' + peptide[2] + 'TAA'
                seq_with_test_pep = current_seq[:end] + \
                    current_seq[end:end+buffer] + \
                    peptide[2] + current_seq[end+buffer:]

                genome_records_dict[gene_header].seq = seq_with_test_pep

                test_data[gene_header] = int(end+buffer)  # add location

        records = [seq_rec for key, seq_rec in genome_records_dict.items()]
        SeqIO.write(records, genome_copy, 'fasta')
        phenotype.genomes[i].genome_file = genome_copy
        phenotype.genomes[i].genome_dict = convert_genome_to_header_dict(
            genome_copy)
        # overwrite the existing genome file with the modified and unmodified

    return test_data


def score_preformance(results_file):
    '''
    Based on the final fasta file output looks for the presence of positive
    test peptides and the absence of negative test peptides. Returns a tuple
    where index one is the percent of positive peptides in final output and
    second index is percent negative peptides in final output. If program
    functioned completly as expected the tuple would = (1, 0).
    '''
    peptide_set = set([str(pep.seq) for pep in SeqIO.parse(results_file, 'fasta')])
    scores = [0, 0]
    with open(TEST_PEPS) as TP:
        reader, totals = csv.reader(TP), [0, 0]
        for row in reader:
            if row[3] == 'P':
                totals[0] += 1
                if row[1] in peptide_set:
                    scores[0] += 1
            else:
                totals[1] += 1
                if row[1] in peptide_set:
                    scores[1] += 1
                    
    return scores[0] / totals[0], scores[1] / totals[1]


def check_gff_file(gff_file, test_data):
    '''
    Checks to make sure that no coding regions overlap with peptide insertions.
    If any overlap is found removes that coding region from the gff file and
    rewrites.
    '''
    good_rows = []
    with open(gff_file) as gff:
        reader = csv.reader(gff, delimiter='\t')
        for row in reader:
            if row[0] in test_data:
                p_s = int(test_data[row[0]])
                if p_s >= int(row[3]) and p_s <= int(row[4]):
                    print('removed region!!!')
                    continue
                else:
                    good_rows.append(row)

    with open(gff_file, 'w') as gff_write:
        writer = csv.writer(gff_write, delimiter='\t')
        for row in good_rows:
            writer.writerow(row)


def clean_up_genome_copies(phenotypes):
    '''
    Removes genomes that have inserted test peptides after the completion
    of a run in testing mode.
    '''
    for p in phenotypes:
        for g in p.genomes:
            if '.copy' in g.genome_file:
                os.remove(g.genome_file)
