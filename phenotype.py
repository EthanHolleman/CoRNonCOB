class Phenotype():
    
    def __init__(self, phenotype):
        self.phenotype = phenotype
        self.genomes = []
        self.conserved_seqs = None
        pass
    
    def add_genome(self, genome):
        '''
        Takes in a genome object and if it has the same phenotype as
        the phenotype instance allows that genotype to be added to the
        phenotypes list of genomes (self.genomes)
        
        :param genome: Genome object
        '''
        if genome.phenotype != self.phenotype:
            return 1
        else:
            self.genomes.append(genome)            
        
    
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

            return shared_seqs, unique_seqs  # want unique to A and to B
            