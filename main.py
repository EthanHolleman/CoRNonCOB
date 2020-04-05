from cornoncob.args_reader import get_args
from cornoncob.phenotype import Phenotype
from cornoncob.io_utils import if_not_exists_make


PROKKA = '/home/ethan/prokka/bin/./prokka'

def main():
    args = get_args()  # parse arguements from the command line
    run_dir = if_not_exists_make(args.o, args.n)
    # make a directory in the output location using the run name
    phenotypes = [Phenotype(pheno_dir, run_dir)
                  for pheno_dir in (args.p1, args.p2)]
    
    for p in phenotypes:
        p.pull_peptides(args.k)
        p.get_conserved_sequences()
    
    phenotypes[0].compare_phenotype(phenotypes[1])
    # phenotype objects are created and directory strucutre is set up


if __name__ == '__main__':
    main()
