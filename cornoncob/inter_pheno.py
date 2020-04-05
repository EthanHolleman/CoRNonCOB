from cornoncob.io_utils import parse_cdhit_output_file

# functions for comparing two phenotypes to get unique peptides

def get_unique_peptides(phenotype_a, phenotype_b, twod_clstr_file):
    '''
    Given two phenptype objects and the cd-hit-2d clstr file produced from
    comparing the phenotypes conserved peptide libraries returns two lists.
    The first being a list of peptides that are uniquely conserved in 
    phenotype_a and the second being one of peptides that are uniqely
    conserved in phenotype_b
    '''
    clstrs = parse_cdhit_output_file(twod_clstr_file)
    
    
    
    
    