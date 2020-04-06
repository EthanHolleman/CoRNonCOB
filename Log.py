from phenotype import Phenotype 

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

    def get_phenotype_parameters(self,phenotype):
        '''
        Takes phenotype object, and writes its parameters to a log file
        '''
        self.write_string('Phenotype: {}\tNum Genomes: {}\tOutput_dir: {}\n'.format(phenotype.phenotype,
len(phenotype.genomes),phenotype.output_dir))

