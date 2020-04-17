from cornoncob.args_reader import get_args
from cornoncob.phenotype import Phenotype
from cornoncob.io_utils import if_not_exists_make
from cornoncob.inter_pheno import get_unique_peptides

from cornoncob import TEST_PEPS  # set in __init__.py
from cornoncob.tests.test_utils import read_test_peps, insert_test_peptides_into_all_phenotypes, score_preformance, check_gff_file
from cornoncob.phenotype import Phenotype
from cornoncob.genome import Genome

from cornoncob.tests.test_utils import check_gff_file

import os
from Bio import SeqIO


def main():

    args = get_args()  # parse arguements from the command line
    run_dir = if_not_exists_make(args.o, args.n)
    # make a directory in the output location using the run name
    phenotypes = [Phenotype(pheno_dir, run_dir)
                  for pheno_dir in (args.p1, args.p2)]

    if args.test:
        test_data = insert_test_peptides_into_all_phenotypes(phenotypes[0], phenotypes[1],
                                                 prokka_exec=args.k)

    for i, p in enumerate(phenotypes):
        for genome in p.genomes:
            genome.make_gene_predictions(path_to_exec=args.k)
            if args.test:
                check_gff_file(genome.gene_prediction_file, test_data[i])
            genome.get_non_coding_regions()
            genome.translate_non_coding_seqs()
            genome.write_peptides_to_fasta_file()
        
        p.get_conserved_sequences()  # compare peptides cross genomes in phemo
        

    unique_seqs = get_unique_peptides(phenotypes[0], phenotypes[1])
    unique_peps_path = os.path.join(run_dir, 'unique_peps.fasta')
    with open(unique_peps_path, 'w') as upp:
        for pep in unique_seqs:
            upp.write(f'>{pep[0]}\n{pep[1]}\n')
    
    if args.test:  # TODO: write to log file
        print(score_preformance(unique_peps_path))
    

    # then if this is a test check to see if we can find the
    # test peptides in the unqiue output

    # phenotype objects are created and directory strucutre is set up

if __name__ == '__main__':
    main()
