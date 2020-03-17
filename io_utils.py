import csv
import os
from Bio import SeqIO


def if_not_exists_make(parent_dir_path, child_dir_name):
    '''
    Checks if a child directory exists within a parent directory. If the child
    directory does not exist makes it. Returns the path to the child directory
    whether it was made or not. Does not search the parent directory
    recursively.

    :param: parent_dir_path: String; Path to the parent directory
    :param: child_dir_name: String; Basename of directory to check if exists
    '''

    child_dir_path = os.path.join(parent_dir_path, child_dir_name)
    if not os.path.exists(child_dir_path):
        os.mkdir(child_dir_path)
    return child_dir_path


def filter_prakka_files(prokka_results_path, *args):
    '''
    Function for returning specific types of files from a directory where
    prokka results have been written. Select the filetypes that should be
    returned by including the file extensions without the . as arguements.

    :params: prokka_results_path: String; Path to directory containing prokka results
    :params: *args: Strings; File extensions for the types of files to be returned
    '''

    args, selected_files = set(args), []
    prokka_results_files = [os.path.join(
        prokka_results_path, f) for f in os.listdir(prokka_results_path)]
    # get all file names in prokka_results_path directory and join with the
    # path to the prokka_results_path to produce list of absolute file paths

    for prokka_file in prokka_results_files:
        if prokka_file.split('.')[-1] in args:
            selected_files.append(prokka_file)
    return selected_files


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


def parse_gff(gff_path, *args, header=False):
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
            if header:  # skip header if exists
                next(reader)
                
            for row in reader:
                if row[0][0] == '#':  # comment in file
                    continue
                else:
                    try:
                        gff_data.append([row[i] for i in args])
                    except IndexError:
                        continue
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