import os
import subprocess

from cornoncob.io_utils import parse_cdhit_output_file, parse_cdhit_output_file

# functions for comparing two phenotypes to get unique peptides


def cd_hit_twod(pep_a, pep_b, output_dir, filename='compare_phenos', s='0.90'):
    '''
    Run cd-hit 2d to compare the peptides in two fasta files. The comparison is
    one way based on the first input so cdhit looks for peptides in phenotype
    b that are similar to phenotype a. Run it with the inputs in reverse to get
    the comparison in the other direction.
    '''
    output_file = os.path.join(output_dir, filename)
    cd_hit_cmd = [
        'cdhit-2d', '-i', pep_a.conserved_seqs, '-i2', pep_b.conserved_seqs, '-o',
        output_file, '-d', '0', '-p', '1', '-s', s
    ]
    subprocess.call(cd_hit_cmd)

    return f'{output_file}.clstr'


def get_unique_peptides(phenotype_a, phenotype_b):
    '''
    Returns peptide headers that are unique to phenotype_a based on the rule
    that a unique peptide will be in a cluster of size one in the cd-hit 2d
    output. This assumption has held so far but probably needs a double check.
    '''
    unique_sequences = []
    clstrs = parse_cdhit_output_file(cd_hit_twod(phenotype_a, phenotype_b, phenotype_a.output_dir))
    for clstr in clstrs:
        if len(clstr) == 1:
            unique_sequences.append((clstr[0][2], phenotype_a.peptide_dict[clstr[0][2]]))
    return unique_sequences
    # this is just headers but actually want headers and sequences need to use the
    # phenotype dictionary to recover that kind of stuff
