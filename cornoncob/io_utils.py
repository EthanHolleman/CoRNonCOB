import csv
import os
from Bio import SeqIO


def if_not_exists_make(parent_dir_path, child_dir_name):
    '''
    Checks if a child directory exists within a parent directory. If the child
    directory does not exist makes it. Returns the path to the child directory
    whether it was made or not. Does not search the parent directory
    recursively.

    :param parent_dir_path: String; Path to the parent directory
    :param child_dir_name: String; Basename of directory to check if exists
    '''

    child_dir_path = os.path.join(parent_dir_path, child_dir_name)
    if not os.path.exists(child_dir_path):
        os.mkdir(child_dir_path)
    return child_dir_path


def filter_prokka_files(prokka_results_path, *args):
    '''
    Function for returning specific types of files from a directory where
    prokka results have been written. Select the filetypes that should be
    returned by including the file extensions without the . as arguements.

    :param prokka_results_path: String; Path to directory containing \
    prokka results
    :param *args: Strings; File extensions for the types of files to be returned
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


def convert_genome_to_header_dict(genome_path, input_format='fasta',
                                  record_format=list):
    '''
    Reads in a file type (likely fasta) and creates a dictionary where keys
    are the sequence headers and the values are lists of nucleotides. Each
    index in the list will hold one nucleotide.
    :param: genome_path: String. Path to the file where genomic \
    sequences are stored
    :param: format: String. File format. Default = fasta
    '''
    records = SeqIO.parse(genome_path, input_format)
    return {record.description: record_format(record.seq) for record in records}


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
    DEPRECIATED

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


def parse_cdhit_record(record, genome_id=True):
    '''
    Parses a single record representing a peptide sequence in a cdhit output
    file. Returns a tuple of that rows contents with extra characters and
    white spaces cleaned up. Tuple contains the index of the sequence in the
    cluster, the length of the sequence, the header of the sequence and the
    alignment information in that order. If genome_id is True meaning a genome
    id has been appended the end of each header this function will add that
    on to the end of the tuple.

    Currently assumes that cd-hit was run with -d 0 so full header is included
    in the output and -p 1 so alignment statistics are printed to the output
    file. 
    '''
    index = record[0]
    length = record.split('\t')[0][:2]
    header, alignment = record.split(',')[1].split('...')
    header = header.split('>')[-1]
    alignment = alignment.strip()
    if alignment != '*':  # is not a representative sequence
        alignment = alignment.split(' ')[-1]

    pretag = [index, length, header, alignment]
    if genome_id:
        tag = header.split('_')[-1].split('...')[0]
        return tuple(pretag + [tag])
    else:
        return tuple(pretag)


def write_unique_seqs(run_dir, unique_seqs, file_name='unique_peps.fasta'):
    '''
    Write peptides from unique_seqs list to fasta output.

    :param: run_dir: String. Path of the current run directory.
    :param: unique_seqs. List. List of unique peptide sequences. Represents \
        the final output of the program.
    '''
    headers = set({})
    unique_peps_path = os.path.join(run_dir, file_name)
    with open(unique_peps_path, 'w') as upp:
        for pep in unique_seqs:
            if pep[0] not in headers:
                upp.write(f'>{pep[0]}\n{pep[1]}\n')
                headers.add(pep[0])
    return unique_peps_path


def parse_cdhit_output_file(clstr_file):
    '''
    Iterates through a cd-hit output file and returns a list of clusters.
    The index of the cluster is equal to the cluster number. Each index contains
    the sequences that cd-hit assigned to that cluster represented as a list
    of tuples. Over data structure looks like [[()]]. Formating of the data
    inside the tuples which represent the actual sequences is determined
    by the parse_cdhit_record function.

    :param clstr_file: String. Path to cd-hit output .clstr file.
    '''
    clusters = []  # store all clusters in list cluster num is index in list
    i = -1
    with open(clstr_file) as clstr:
        while clstr:
            cur_line = clstr.readline().strip()
            if cur_line:
                if cur_line[0] == '>':  # new cluster
                    clusters.append([])
                    i += 1
                else:
                    clusters[i].append(parse_cdhit_record(cur_line))
            else:
                break
        return clusters
