from cornoncob.phenotype import Phenotype
from cornoncob.genome import Genome 

class Log():
    '''
    Methods for writing to log file
    
    '''
    def __init__(self, logfile):
        self.logfile = open(logfile, 'w')
        
    def write_string(self, string):
        self.logfile.write('{}\n'.format(string))
        if args.test:
        # test I use to see if program is in testing mode
            self.logfile.write('Progam is in testing mode')
        
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
        self.write_string('\nGenome: {}\tPhenotype: {}\tNumber of non-coding sequences: {}\n'.format(genome.genome_id,
                                                                                                     genome.phenotype, len(genome.non_coding_sequences)))
        

