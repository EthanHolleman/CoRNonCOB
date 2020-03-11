'''
Put helper functions for reading and parsing files in this file.
'''
import csv

def parse_gff(gff_path, *args):
    '''
    Given a path to a gff file returns as a list of lists the information
    contained in the field indicies given via *args. All args should be
    integers.
    
    1 	sequence  # subtract 1 from these!!!!!
    2 	source
    3 	feature
    4 	start
    5 	end
    6 	score
    7 	strand
    8 	phase
    9 	attributes   
    '''   
    gff_data = []
    if args:
        with open(gff_path) as gff:
            reader = csv.reader(gff, delimiter='\t')
            for row in reader:
                gff_data.append([row[i] for i in args])
    return gff_data


def genome_parser(genome_path):
    '''
    Reads a genome file and returns the genome as a list where each index
    is one nucleotide. list[0] = first nucleotide
    '''
    pass