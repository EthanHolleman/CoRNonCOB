from cornoncob.phenotype import Phenotype
from cornoncob.genome import Genome
import time 

class Log():
    '''
    Methods for writing to log file
    
    '''
    def __init__(self, logfile):
        self.logfile = open(logfile, 'w')
        self.start_time = time.time()
        
    def write_string(self, string):
        self.logfile.write('{}'.format(string))
        
    def close_log(self):
        self.logfile.close()

    def get_phenotype_parameters(self,phenotype):
        
        '''

        Method that takes phenotype object,
        and writes its parameters to a log file
        
        '''
        
        self.write_string('Phenotype: {}\tNum Genomes: {}\tOutput_dir: {}\n'.format(phenotype.phenotype,
len(phenotype.genomes),phenotype.output_dir))

    def get_genotype_parameters(self, genome):
        
        '''

        Method that takes in genome object,
        and writes number of nc regions to file.
        
        '''
        self.write_string('Genome: {}\tPhenotype: {}\tNumber of non-coding sequences: {}'.format(genome.genome_id,
                                                                                                     genome.phenotype, len(genome.non_coding_sequences)))
    def get_number_conserved_peptides(self,phenotype):

        '''

        Method that takes phenotype object,
        and writes number of conserved peptides
        
        '''
        cp = 0
        f = open(phenotype.conserved_seqs, 'r')
        for line in f:
            if line[0] == '>':
                cp = cp + 1

        
        self.write_string('\nThere are {} conserved peptides in phenotype {}'.format(cp, phenotype.phenotype))

    def get_number_unique_peptides(self, directory):

        '''

        Method that takes the directory of the fasta file,
        and writes number of unique peptides
        
        '''
        up = 0
        f2 = open(directory, 'r')
        for line in f2:
            if line[0] == '>':
                up = up + 1

        self.write_string('\nThere are {} unique peptides'.format(up))

    def calculate_execution_time(self):
        '''

        Method that records execution time of program
         
        '''

        elapsed_time = time.time()- self.start_time

        self.write_string('\nExecution time: %s seconds' % elapsed_time)

