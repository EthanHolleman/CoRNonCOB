from cornoncob.args_reader import get_args
from cornoncob.phenotype import Phenotype
from cornoncob.io_utils import if_not_exists_make
from cornoncob.inter_pheno import get_unique_peptides

from cornoncob import TEST_PEPS  # set in __init__.py
from cornoncob.tests.test_utils import read_test_peps, run_baseline_prokka, insert_test_peps
from cornoncob.phenotype import Phenotype
from cornoncob.genome import Genome

from cornoncob.tests.test_utils import check_gff_file

import os


def main():

    args = get_args()  # parse arguements from the command line
    run_dir = if_not_exists_make(args.o, args.n)
    # make a directory in the output location using the run name
    phenotypes = [Phenotype(pheno_dir, run_dir)
                  for pheno_dir in (args.p1, args.p2)]

    if args.test:
        test_data = insert_test_peps(
            phenotypes[0], TEST_PEPS, seed=100, prokka_exec=args.k)

    for p in phenotypes:
        p.pull_peptides(prokka_exec=args.k)  # convert noncoding -> peptides
        p.get_conserved_sequences()  # compare peptides cross genomes in phemo

    unique_seqs = get_unique_peptides(phenotypes[0], phenotypes[1])
    with open(os.path.join(run_dir, 'unique_peps.fasta'), 'w') as unique_peps:
        for pep in unique_seqs:
            unique_peps.write(f'>{pep[0]}\n{pep[1]}\n')

    # then if this is a test check to see if we can find the
    # test peptides in the unqiue output

    # phenotype objects are created and directory strucutre is set up

if __name__ == '__main__':
    main()
