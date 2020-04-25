from cornoncob.io_utils import parse_gff, if_not_exists_make, convert_genome_to_header_dict
from cornoncob.phenotype import Phenotype
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import csv
import sys
import subprocess
import os

from cornoncob import TEST_PEPS, TEST_KILLERS, TEST_NICE
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


def make_test_phenotypes(run_dir):
    # just for testing make based on test data
    pos_test_peps, neg_test_peps = read_test_peps(TEST_PEPS)
    pos_test_genomes, pos_in_sites = insert_test_peptides_into_phenotype(TEST_KILLERS,
                                                                        pos_test_peps)
    neg_test_genomes, neg_in_sites = insert_test_peptides_into_phenotype(TEST_NICE,
                                                                         neg_test_peps)
    return [Phenotype(pos_test_genomes, run_dir),
            Phenotype(neg_test_genomes, run_dir)], pos_in_sites, neg_in_sites
    
    # return list of the phenotypes set up based on the inserted test
    # peptide files 
    # now just need to after run prokka on these in test mode run through and
    # remove the coding regions which contain the insterted test peptides
    # by modifying the gff file before getting non-coding regions
    


def insert_test_peptides_into_phenotype(pheno_dir, test_peptide_list, stop=True,
                                        seed=100, prokka_exec='prokka'):
    '''
    Inserts the given test_peptides from the test_peptide_list into the
    given phenotype at random contigs. Returns a dictionary of contig headers
    and locations of peptide insertions.
    '''
    np.random.seed(seed)
    test_data = {}
    num_test_peps = len(test_peptide_list)
    insertion_sites_dict = {}  # filepath to list of instertion sites
    test_dir = if_not_exists_make(pheno_dir, os.path.basename(pheno_dir) +'_test_genomes')
    
    for genome_file in os.listdir(pheno_dir):
        insertion_sites = {}
        abs_genome_file = os.path.join(pheno_dir, genome_file)
        if os.path.isdir(abs_genome_file):
            continue
        print(abs_genome_file)
        genome_records_dict = SeqIO.to_dict(
                SeqIO.parse(abs_genome_file, 'fasta')
            )  # turn genome file into dictionary by header
        random_record_headers = np.random.choice(
                len(genome_records_dict.keys()), num_test_peps, replace=False
            )
        random_record_headers = [list(genome_records_dict.keys())[i] for i in random_record_headers]
        # pick random headers to insert into
        for j, random_header in enumerate(random_record_headers):
            peptide = 'TAA' + test_peptide_list[j][2] + 'TAA'
            seq = str(genome_records_dict[random_header].seq)
            insertion_start = np.random.randint(20, len(seq))
            insertion_end = insertion_start + len(peptide)
            inserted_seq = seq[:insertion_start] + peptide + seq[insertion_start:]
            
            genome_records_dict[random_header].seq = Seq(inserted_seq)
            # reassing the seq with inserted peptide to the header
            insertion_sites[random_header] = (insertion_start,
                                    insertion_end)
        
        genome_copy = f'{genome_file}.copy'  # genome copy file name
        genome_copy_file_path = os.path.join(test_dir, genome_copy)  # full path
        print(genome_copy_file_path)
        records = [seq_rec for key, seq_rec in genome_records_dict.items()]
         # convert dict back to just records
        SeqIO.write(records, genome_copy_file_path, 'fasta')
        print('wrote dat file', len(records))
        # write the copy file
        insertion_sites_dict[genome_copy_file_path] = insertion_sites
    
    return test_dir, insertion_sites_dict
    # now need to make the phenotypes using this new test data


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
                    print(row[1])
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
    
    Test data is now a dictionary with the header as key and the insertion
    start and end sites as values
    test data should just be the dictory of headers
    '''
    good_rows = []
    with open(gff_file) as gff:
        reader = csv.reader(gff, delimiter='\t')
        for row in reader:
            if row[0] in test_data:
                gene_start, gene_end = int(row[3]), int(row[4])
                p_s, p_e = test_data[row[0]]  # start position end position
                if (p_s < gene_start and p_e > gene_start) or (p_s < gene_end and p_e > gene_end) or (p_s > gene_start and p_e < gene_end):
                    print('removed region', p_s, p_e, gene_start, gene_end)
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
