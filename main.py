import os
from Bio import SeqIO

from cornoncob import TEST_PEPS  # set in __init__.py
from cornoncob.args_reader import get_args
from cornoncob.genome import Genome
from cornoncob.inter_pheno import get_unique_peptides
from cornoncob.io_utils import if_not_exists_make, write_unique_seqs
from cornoncob.phenotype import Phenotype
from cornoncob.tests.test_utils import (
    check_gff_file, insert_test_peptides_into_all_phenotypes, read_test_peps,
    score_preformance, clean_up_genome_copies)
from cornoncob.chemical_props import write_peptide_properties
from cornoncob.Log import *


def main():

    args = get_args()  # parse arguements from the command line
    run_dir = if_not_exists_make(args.o, args.n)
    log = Log(run_dir + '.log')
    # make a directory in the output location using the run name
    phenotypes = [Phenotype(pheno_dir, run_dir)
                  for pheno_dir in (args.p1, args.p2)]
    for phenotype in phenotypes:
        log.get_phenotype_parameters(phenotype)
    

    if args.test:  # program is run in test mode so run initial prokka to insert test peps
        log.write_string('CoRNonCOB Testing Mode')
        test_data = insert_test_peptides_into_all_phenotypes(phenotypes[0], phenotypes[1],
                                                 prokka_exec=args.k)

    for i, p in enumerate(phenotypes):
        for genome in p.genomes:
            genome.make_gene_predictions(path_to_exec=args.k)
            if args.test:  # make sure no test peps in coding regions if test 
                check_gff_file(genome.gene_prediction_file, test_data[i])
            genome.get_non_coding_regions()
            genome.translate_non_coding_seqs()
            genome.write_peptides_to_fasta_file()
        
        p.get_conserved_sequences()  # compare peptides cross genomes in phemo
        log.get_number_conserved_peptides(p)
        

    unique_seqs = get_unique_peptides(phenotypes[0], phenotypes[1])
    unique_peps_path = write_unique_seqs(run_dir, unique_seqs)
    chemical_props = write_peptide_properties(unique_peps_path, run_dir)
    log.get_number_unique_peptides(unique_peps_path)
    
    if args.test:  # TODO: write to log file
        clean_up_genome_copies(phenotypes)
        log.get_number_unique_peptides(unique_peps_path)
        log.get_score_preformance(score_preformance(unique_seqs_path))
        
    log.calculate_execution_time()
    
    # then if this is a test check to see if we can find the
    # test peptides in the unqiue output

    # phenotype objects are created and directory strucutre is set up

if __name__ == '__main__':
    main()
