import subprocess
import os
import Bio

from cornoncob.io_utils import *
from cornoncob.sequence import NoncodingSeq


class Genome():
    '''
    Representation of an individual bacterial genome. Stores its NoncodingSeq files
    and has methods for predicting coding regions and then retreiving
    non-coding regions

    :param genome_file: string, filepath to fasta formated file of \
        the complete bacterial genome
    :param phenotype: string, string that acts as a key for identifying \
        the bacterial phenotype
    :param gene_prediction_file: string, filepath to gff generated via prokka.\
        Default = None
    :param non_coding_file: string, filepath to fasta of non-coding regions.\
        Default = None
    '''

    def __init__(self, genome_file, phenotype, pheno_dir,
                 gene_prediction_file=None, non_coding_file=None):

        self.genome_file = genome_file
        self.genome_id = str(hash(genome_file))[:5]
        self.phenotype = phenotype
        self.gene_prediction_file = gene_prediction_file
        self.non_coding_file = non_coding_file
        self.non_coding_seqs = []
        self.output_dir = if_not_exists_make(
            pheno_dir, os.path.basename(genome_file))
        self.genome_dict = None
        # potentially change to dictionary to allow non coding seq
        # look up by start location if that is needed later

    def __repr__(self):
        '''
        Determines what is returned when a genome object is printed.
        '''
        s = ['Phenotype:', str(self.phenotype), 'Genome File:', str(self.genome_file),
             'Gene Prediction File:', str(self.gene_prediction_file),
             'Non-coding File:', str(self.non_coding_file),
             'Num non-coding seqs:', str(len(self.non_coding_seqs))]
        return ' '.join(s)

    def __len__(self):
        return len(self.non_coding_seqs)

    def __getitem__(self, peptide_header=None):
        if peptide_header:
            pass

    def make_gene_predictions(self, threads=2, path_to_exec='prokka',
                              results_dir_name='prokka_results'):
        '''
        This function will run prokka to find predicted gene locations for
        the instance of Genome it is called on.
        '''
        # TODO: Get this function woring with correct variable names
        # and output paths organized
        # will need to make new folder for prokka results for each genome
        # that runs gene prediction
        prokka_dir = if_not_exists_make(self.output_dir, results_dir_name)

        input_file = self.genome_file
        

        cmd = [path_to_exec, '--outdir', prokka_dir,
               '--cpus', str(threads), '--force', '--fast', input_file]
        subprocess.call(cmd)
        output_file = filter_prokka_files(
            self.output_dir + '/prokka_results', 'gff')[0]
        self.gene_prediction_file = output_file
    # run gene prediction methods here

    def get_non_coding_regions(self):
        # may have to modify if genomes are not all on one chromosome
        '''
        Using the genome file associated with an instance of a specific genome
        object and the gene prediction gff file to extract the non-coding
        regions in the genome. Will write the non-coding strings as a fasta
        file and store the path to that file in the non_coding_file
        attribute.
        '''
        self.genome_dict = convert_genome_to_header_dict(self.genome_file)
        coding_regions = parse_gff(self.gene_prediction_file, 0, 3, 4)
        # read the gff file for predicted coding regions

        coding_regions_dict = {}

        # make dictionary based on NoncodingSeq headers. Each header key is mapped
        # to a list of tuples. Each tuple contains the start and stop location
        # for a predicted gene

        for header, start, stop in coding_regions:
            start, stop = int(start), int(stop)
            if header in coding_regions_dict:
                coding_regions_dict[header].append((start, stop))
            else:
                coding_regions_dict[header] = [(start, stop)]

        # sort the coding regions for eac key (header) in coding_regions_dict
        for header, cordinates in coding_regions_dict.items():
            coding_regions_dict[header] = sorted(
                cordinates, key=lambda x: x[0])
            # sort by start position of the predicted coding region

        for header, coding_regions_list in coding_regions_dict.items():
            # iterate through all pairs of headers and coding regions
            # get the sequence for the current header
            
            working_seq = self.genome_dict[header]
            cur_non_coding_string = ''
            i, j = 0, 0
            start, stop = coding_regions_list[j]
            while True:
                
                if i < start:  # at position before gene
                    cur_non_coding_string += working_seq[i]
                    i += 1
                else:  # in a gene
                    description = '{}.{}'.format(header, i)
                    self.non_coding_seqs.append(
                        NoncodingSeq(description, cur_non_coding_string, i))
                    cur_non_coding_string = ''
                    i = stop  # jump to end of the gene and get cordinates
                    # for next coding region

                    j += 1
                    if j < len(coding_regions_list):
                        prev_stop = stop
                        start, stop = coding_regions_list[j]
                        if start < prev_stop and j+1 < len(coding_regions_list):
                            j += 1  # skip to next since overlapping region found
                            start, stop =coding_regions_list[j]      
                    else:
                        self.non_coding_seqs.append(
                            NoncodingSeq(description, ''.join(working_seq[i:]), i))
                        break
        not_in_coding = 0
        for fasta_header in self.genome_dict:  # if no coding regions in record
            if fasta_header not in coding_regions_dict:
                not_in_coding += 1
                self.non_coding_seqs.append(
                    NoncodingSeq(fasta_header,
                                ''.join(self.genome_dict[fasta_header]),
                                0))  # add the whole record
        print(not_in_coding, ' entries ', 'of ', len(self.genome_dict.keys()))

    def translate_non_coding_seqs(self):
        '''
        Iterates through the non coding NoncodingSeqs stored in non_coding_seqs
        and calls translation_siz_shooter to translate the nucleotide
        NoncodingSeq into all six reading frames.
        '''
        for ncs in self.non_coding_seqs:
            try:
                ncs.translation_six_shooter()
            except Bio.Data.CodonTable.TranslationError as e:
                continue

    def write_peptides_to_fasta_file(self, output_dir=None,
                                     stop_codon_symbol='*', min_len=6):
        '''
        After the noncoding peptides have been translated this function is
        called to write those peptides to a file so a motif finding or
        clustering program can find conserved and or differential NoncodingSeqs.
        Currently the headers for each peptide are pretty minimal and will be
        reworked to include more info soon.

        :param output_dir: String. Path to directory where fasta file will be written
        :param stop_codon_symbol: Char. Symbol used for stop codon in translated seqs. Default = '*'
        :param min_len: Int. Minimum length required for peptide to be written.
        '''
        if output_dir == None:
            output_dir = self.output_dir

        file_basename = '{}_peptides'.format(
            os.path.basename(self.genome_file))
        self.non_coding_file = os.path.join(output_dir, file_basename)

        with open(self.non_coding_file, 'w') as fp:
            for noncoding_seq in self.non_coding_seqs:
                # get noncoding seq objects
                for i, frame in enumerate(noncoding_seq.frames):
                    # get individual reading fram
                    peptides = str(frame).split(stop_codon_symbol)
                    for j, pep in enumerate(peptides):
                        if len(pep) >= min_len:
                            fp.write('>{}_{}_{}_{}\n'.format(
                                noncoding_seq.description, i, j, self.genome_id))
                            # i is the frame number and j in the peptide index in that frame
                            fp.write(pep + '\n')
                    # write the individual peptides

    # read in the genome file
    # get positions of coding regions
    # extract the non-coding regions and write to a file
