from modlamp.descriptors import PeptideDescriptor, GlobalDescriptor
from Bio import SeqIO
import csv
import os


def get_sequence_dict(fasta_path):
    '''
    Turns a fasta file into a dictionary where the keys are the headers and the
    values are SeqRecord objects.
    '''
    return SeqIO.to_dict(SeqIO.parse(fasta_path, 'fasta'))


def calculate_peptide_props(fasta_dict):
    '''
    Give a sequence_dictionary (made from get_sequence_dict) returns a
    list of dictionaries. Each dictionary has type of chemical property as the
    keys and the calculated value for that property as the value. Designed to
    be written to a csv file using DictWriter.
    '''
    property_list = []
    for header in fasta_dict:
        s = str(fasta_dict[header].seq)
        t = GlobalDescriptor([s])
        t.calculate_all()
        d = dict(zip(t.featurenames, t.descriptor[0]))
        d['Peptide_name'] = header
        property_list.append(d)
    return property_list


def write_peptide_properties(input_fasta, run_dir,
                             file_name='chemical_props.csv'):
    '''
    Given a fasta file with peptides and the current run dir writes as csv
    file with calculated chemical properties for each peptide. Each row in the
    csv output is one peptide from the input_fasta file.
    '''
    output_csv = os.path.join(run_dir, file_name)
    property_dict = calculate_peptide_props(get_sequence_dict(input_fasta))
    with open(output_csv, 'w') as output:
        writer = csv.DictWriter(output, fieldnames=property_dict[0].keys())
        writer.writeheader()
        for peptide in property_dict:
            writer.writerow(peptide)
    return output_csv
