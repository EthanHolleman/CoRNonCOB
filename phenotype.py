class Phenotype():
    
    def __init__(self, phenotype):
        self.phenotype = phenotype
        self.genomes = []
        self.conserved_seqs = None
        pass
    '''
    Potential class for holding the genomes that are under each
    phenotype. Could be useful for organizing functions that compare
    specific regions of genomes between phenotypes.
    
    Or it might be better to have a class that just holds all the genomes
    in one data structure and can retrieve them based on phenotype that way
    all comparisons can be done with methods just in that one object. 
    
    Or this set up could have functions that take in another instance of a
    phenotype and do the comparisons between the instance the method is called
    on and the phenotype object that is passed in.
    
    Liking last option the best so far
    
    '''
    
    def add_genome(self, genome):
        '''
        Takes in a genome object and if it has the same phenotype as
        the phenotype instance allows that genotype to be added to the
        phenotypes list of genomes (self.genomes)
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
            