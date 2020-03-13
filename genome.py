from io_utils import parse_gff
from io_utils import convert_genome_to_list
from sequence import Sequence


class Genome():
    '''
    Representation of an individual bacterial genome. Stores its sequence files
    and has methods for predicting coding regions and then retreiving
    non-coding regions

    :param genome_file: string, filepath to fasta formated file of the complete bacterial genome
    :param phenotype: string, string that acts as a key for identifying the bacterial phenotype
    :param gene_prediction_file: string, filepath to gff generated via prokka. Default = None
    :param non_coding_file: string, filepath to fasta of non-coding regions. Default = None
    '''

    def __init__(self, genome_file, phenotype, gene_prediction_file=None, non_coding_file=None):
        self.genome_file = genome_file
        self.phenotype = phenotype
        self.gene_prediction_file = gene_prediction_file
        self.non_coding_file = non_coding_file
        self.non_coding_seqs = []
        # potentially change to dictionary to allow non coding seq
        # look up by start location if that is needed later

    def make_gene_predictions(self):
        pass
    # run gene prediction methods here

    def get_non_coding_regions(self):
        '''
        Using the genome file associated with an instance of a specific genome
        object and the gene prediction gff file to extract the non-coding
        regions in the genome. Will write the non-coding strings as a fasta
        file and store the path to that file in the non_coding_file
        attribute.
        '''
        start_stop_list = parse_gff(self.gene_prediction_file, 3, 4)
        genome = convert_genome_to_list(self.genome_file)
        cur_non_coding_string = ''
        i, j = 0, 0
        start, stop = start_stop_list[j]

        while True:
            if i < start:
                cur_non_coding_string += genome[i]
                i += 1
            else:
                self.non_coding_seqs.append(Sequence(cur_non_coding_string, i)
                cur_non_coding_string = ''
                i = stop
                # need to look into this more because the positions get from
                # the actaul gff file will be base one while i is refering
                # to positions in the genome in base 0
                j += 1
                if j < len(start_stop_list):
                    start, stop = start_stop_list[j]
                else:
                    self.non_coding_seqs.append(''.join(Sequence(genome[i:], i))
                    break

    # read in the genome file
    # get positions of coding regions
    # extract the non-coding regions and write to a file