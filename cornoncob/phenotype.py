import os
import subprocess
import sys

from cornoncob.genome import Genome
from cornoncob.io_utils import (convert_genome_to_header_dict,
                                if_not_exists_make, parse_cdhit_output_file)


class Phenotype():
    '''
    Designed to represent a specific
    bacterial phenotype that may contain many individual genomes of that
    phenotype. Phenotype objects are currently the outermost layer of the
    program.

    :param genome_dir: String. Path to directory containing all \
        genomes that are to be contained in the phenotype instance.
    :param run_dir: String. Path to outermost directory and is supplied \
        by the user. Creating an instance of phenotype will create a new \
        directory within the run_dir.
    :param phenotype: String. Description of the phenotype. If none is \
        provided assumes that the basename of the genome_dir is the \
        phenotype name. 
    '''

    def __init__(self, genome_dir, run_dir, phenotype=None):
        self.genome_dir = genome_dir
        self.phenotype = self.set_phenotype(phenotype)
        self.output_dir = if_not_exists_make(run_dir, self.phenotype)
        self.genomes = self.add_genomes_from_dir(genome_dir, phenotype)
        # store the fasta files as Genome objects in self.genomes
        self.conserved_seqs = None
        self.peptide_dict = {}

    def __repr__(self):
        return 'Phenotype: {}\nNum Genomes: {}\nOutput_dir: {}'.format(self.phenotype,
                                                                       len(self.genomes),
                                                                       self.output_dir)

    def __getitem__(self, id_header_tuple):
        genome_id, seq_header = id_header_tuple
        for genome in self.genomes:
            if genome.genome_id == genome_id:
                if seq_header in genome.genome_dict:
                    return genome.genome_dict[seq_header]
                else:
                    print('Header not found in genome dict', seq_header)
            else:
                print('genome id not found', genome_id)

    def set_phenotype(self, phenotype):
        if phenotype:
            return phenotype
        else:
            return os.path.basename(self.genome_dir)

    def add_genomes_from_dir(self, genome_dir, pheno=None):
        '''
        Takes a directory containing genome files of one phenotype and adds those
        files as genome objects to this phenotype instance's genomes list.

        param: genome_dir: String; path to the genome directory
        param: pheno: String; If not None, use this string as the phenotype 
        '''
        genomes = []
        # deciede what will be used for the phenotype name
        if pheno == None:
            self.phenotype = os.path.basename(genome_dir)
        else:
            self.phenotype = pheno

        if genome_dir:  # not none
            genomes_paths = [os.path.join(genome_dir, genome)
                             for genome in os.listdir(genome_dir)]
            for genome_path in genomes_paths:
                genomes.append(Genome(genome_path, pheno, self.output_dir))

        return genomes

    def pull_peptides(self, prokka_exec='prokka'):
        '''
        Wrapper around methods in the Genome class. This methos iterates
        through all genomes in a phenotype and uses prokka to make gene
        predictions, identifies non-coding regions and then translates
        those regions into six reading frames.
        '''
        for genome in self.genomes:
            genome.make_gene_predictions(path_to_exec=prokka_exec)
            genome.get_non_coding_regions()
            genome.translate_non_coding_seqs()
            genome.write_peptides_to_fasta_file()

    def get_conserved_sequences(self, cdhit_exec='cdhit', s='0.90', con=0.75):
        '''
        :param cdhit_exec: String. Path to cdhit executable default = cdhit
        :param s: String. Num 0-1 sets min length difference between rep seq \
        and subject seqs
        :param con: Int. Percentage of genomes that must have a peptide in a \
        cluster for that cluster to be considered conserved.
        '''
        s = str(s)  # cast in case an int or float gets passed in 

        # get all file paths together
        phenotype_peptides = os.path.join(
            self.output_dir, f'{self.phenotype}_peptides.fasta')
        clstr_file = os.path.join(
            self.output_dir, f'{self.phenotype}_peptides')
        conserved_seq_fasta = os.path.join(
            self.output_dir, f'{self.phenotype}_conserved_peptides.fasta')

        genome_peptides = [
            f'"{genome.non_coding_file}"' for genome in self.genomes]
        # get all of the peptide file paths in one list

        cat_cmd = ['cat'] + genome_peptides + ['>', phenotype_peptides]
        cd_hit_cmd = [cdhit_exec, '-i',
                      phenotype_peptides, '-o', clstr_file, '-d', '0',
                      '-p', '1', '-s', s]

        # concat all individual peptide files
        cat_call = subprocess.call(' '.join(cat_cmd), shell=True)
        cdhit_call = subprocess.call(cd_hit_cmd)  # run cd-hit on cated file

        clusters = parse_cdhit_output_file(clstr_file + '.clstr')
        # cdhit always adds .clstr suffix to results file

        self.peptide_dict = convert_genome_to_header_dict(
            phenotype_peptides, record_format=str)

        conserved_records = []

        for cluster in clusters:
            participating_genomes, rep_seq = set([]), None
            for record in cluster:
                if record[-2] == '*':
                    rep_seq = record
                participating_genomes.add(record[-1])  # add genome id
            if len(participating_genomes) / len(self.genomes) >= con:
                # percentage of genomes that participate in this cluster
                # is greater than or equal to con threshold
                conserved_records.append((rep_seq, len(cluster)))

        with open(conserved_seq_fasta, 'w') as csf:
            for record in conserved_records:  # write the conserved records
                header = record[0][2]
                seq = self.peptide_dict[header]
                csf.write(f'>{header}\n{seq}\n')

        self.conserved_seqs = conserved_seq_fasta

    def assign_headers_to_seqs(self, header_list):
        '''
        Creates a list of tuples with first item is the header of a fasta record
        and the second item is the sequence for that record
        '''
        return [(header, self.peptide_dict[header]) for header in
                header_list if header in self.peptide_dict]
