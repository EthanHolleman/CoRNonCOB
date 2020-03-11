

class Genome():
    
    def __init__(self, genome_file, phenotype):
        self.genome_file = genome_file
        self.phenotype = phenotype
        self.gene_prediction_file = None
        self.non_coding_file = None
    
    
    def make_gene_predictions(self):
        pass
    # run gene prediction methods here
    
    def get_non_coding_regions(self):
        pass
    # read in the genome file
    # get positions of coding regions
    # extract the non-coding regions and write to a file

    def translate_non_coding(self):  # oxymoron?
        pass
        # convert non_coding sequnces into all 6 reading frames
        # not sure if we want to keep track of this if we do becuase stuff
        # like the position would be important then should make a seperate
        # object for sequence or non coding or something like that that
        # can then be translated
    

'''
testing get coding regions. This function works but will need to be adapted
with parsers to pull in the info it needs from actual files

Currently genome is a list with one nucleotide at each index and 
start_stop list is a list of tuples with the start location and end
locations of predicted genes.



def get_non_coding_regions(genome, start_stop_list):
    start, stop = start_stop_list[0]
    non_coding_strings = []
    cur_non_coding_string = ''
    i, j = 0, 0
    start, stop = start_stop_list[j]
    
    while True:
        if i < start:
            cur_non_coding_string += genome[i]
            i += 1    
        else:
            non_coding_strings.append(cur_non_coding_string)
            cur_non_coding_string = ''
            i = stop + 1
            j += 1
            if j < len(start_stop_list):
                start, stop = start_stop_list[j]
            else:
                non_coding_strings.append(''.join(genome[i:]))
                break

    return non_coding_strings
'''


