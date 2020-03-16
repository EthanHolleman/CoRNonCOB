import csv
from Bio import SeqIO



def read_phenotype_map(map_path):
    '''
    If the user has genomes all over the place and does not
    want to aggregate them they can provide a csv file with the
    headers below to create Phenotype objects in the main
    function.

    PHENOTYPE, FILEPATH
    '''
    pass
    # TODO: Write this function and add method in phenotype or main
    # to make phenotype object from this kind of input

    # Should return a dictionary key = filepath, value = phenotype


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


def convert_genome_to_list(genome_path, format='fasta'):
    '''
    Takes in the path to a genome in the format specified by the format
    variable. Genome file is read using SeqIO from Biopython and if the file
    contains multible records the nucleotide sequences are appended to one
    string. String is then returned as a list with each index containing one
    nucleotide.
    '''
    genome_string = ''
    genome_file = SeqIO.parse(genome_path, format)
    for record in genome_file:
        genome_string += record.seq
    return [nuc for nuc in genome_string]
