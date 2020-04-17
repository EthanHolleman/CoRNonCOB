from cornoncob.io_utils import parse_gff, if_not_exists_make, convert_genome_to_header_dict
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import csv
import sys
import subprocess

from cornoncob import TEST_PEPS
sys.path.append(".")


def read_test_peps(TEST_PEPS):
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
    scores = [0, 0]
    p, n = read_test_peps(TEST_PEPS)
    all_peps = p + n
    test_pep_dict = {pep[1]: pep[3] for pep in all_peps}
    records = SeqIO.parse(results_file, 'fasta')
    for r in records:
        str_seq = str(r.seq)
        if str_seq in test_pep_dict:
            if test_pep_dict[str_seq] == 'P':
                scores[0] += 1
            elif test_pep_dict[str_seq] == 'C':
                scores[1] += 1
    num_pos_peps = sum([1 for s, t in test_pep_dict.items() if t == 'P'])

    return scores[0] / num_pos_peps, scores[1] / len(n)
    # returns the percent of positive peptides recovered 1 means got all
    # peptides we expected where returned in the final output
    # and the percent of negative peptides found in the final output. If
    # preformed perfectly we would exptect to this value to be 0 since non
    # of the negative peptides should be found in the final output


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
