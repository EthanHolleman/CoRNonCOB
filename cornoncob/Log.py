from cornoncob.phenotype import Phenotype 

class Log():
    '''
    Methods for writing to log file
    
    '''
    def __init__(self, logfile):
        self.logfile = open(logfile, 'w')
        
    def write_string(self, string):
        self.logfile.write('{}\n'.format(string))
        
    def close_log(self):
        self.logfile.close()

    def get_phenotype_parameters(self, phenotype):
        self.write_string('Phenotype: ' + phenotype.phenotype)
        self.write_string('Genome Directory: ' + phenotype.genome_dir)
        self.write_string('Output Directory: ' + phenotype.output_dir)
        self.write_string('Genomes: ' + phenotype.genomes)
