import os

from genome import Genome


class Phenotype():

    def __init__(self, genome_dir, phenotype=None):
        self.phenotype = phenotype
        self.genome_dir = genome_dir
        self.genomes = []
        self.conserved_seqs = None

        self.add_genomes_from_dir(self.genome_dir, pheno=self.phenotype)

    def add_genome(self, genome_path, phenotype):
        '''
        Given a path to a genome file (fasta file) and a phenotype
        adds an instance of a genome object to the phenotype's
        genomes list attribute.
        '''

        self.genomes.append(Genome(genome_path, phenotype))

    def add_genomes_from_dir(self, genome_dir, pheno=None):
        '''
        Takes a directory containing genome files of one phenotype and adds those
        files as genome objects to this phenotype instance's genomes list.

        param: genome_dir: String; path to the genome directory
        param: pheno: String; If not None, use this string as the phenotype 
        '''

        # deciede what will be used for the phenotype name
        if pheno == None:
            self.phenotype = os.path.basename(genome_dir)
        else:
            self.phenotype = pheno

        if genome_dir:  # not none
            genomes_paths = [os.path.join(genome_dir, genome)
                             for genome in os.listdir(genome_dir)]
            for genome_path in genomes_paths:
                self.add_genome(genome_path, self.phenotype)

    def get_conserved_sequences(self):
        '''
        Using some metrics and methods pulls out the conserved non-coding
        sequences from a collection of genomes stored in this phenotype
        instance.
        '''
        pass

    def compare_conserved_sequences(self, other_phenotype):
        if self.phenotype != other_phenotype:
            # do the comparison since same phenotype
            # will need to do some kind of statistical test here
            # good to provide p values and confidence intervals

            return shared_seqs, unique_seqs  # want unique to A and to B
