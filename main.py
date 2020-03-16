from args import get_args
from phenotype import Phenotype

def main():
    run_args = get_args()  # parse arguements from the command line


    #TODO: Add someway of specifying phenotype name if not given by dir names
    pheno_a = Phenotype(run_args.p1)  # set up Phenotypes
    pheno_b = Phenotype(run_args.p2)

    # if phenotypes are not the dir they are in possible
    # have user supply csv file that acts as table
    # for what the phenotypes are and where that genome can
    # be located





    # maybe is have time add n number of phenotype comparisons
    # using ANOVA methods


    pass

if __name__ == '__main__':
    main()