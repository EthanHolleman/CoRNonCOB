import os
import subprocess

from genome import Genome
from io_utils import if_not_exists_make, parse_cdhit_output_file


# TODO
# phenotype objects should create their own directory within the given
# run directory then potentially hand that directory off to genome objects
# when they are created and in genome __init__ functions genome instance
# will make its own directory with the phenotype dir as the parent directory

class Phenotype():
    '''
    Designed to represent a specific
    bacterial phenotype that may contain many individual genomes of that
    phenotype. Phenotype objects are currently the outermost layer of the
    program.

    :param genome_dir: String. Path to directory containing all \
        genomes that are to be contained in the phenotype instance.
    :param run_dir: String. Path to outermost directory and is supplied \
        by the user. Creating an instance of phenotype will create a new \
        directory within the run_dir.
    :param phenotype: String. Description of the phenotype. If none is \
        provided assumes that the basename of the genome_dir is the \
        phenotype name. 
    '''

    def __init__(self, genome_dir, run_dir, phenotype=None):
        self.genome_dir = genome_dir
        self.phenotype = self.set_phenotype(phenotype)
        self.output_dir = if_not_exists_make(run_dir, self.phenotype)
        self.genomes = self.add_genomes_from_dir(genome_dir, phenotype)
        # store the fasta files as Genome objects in self.genomes
        self.conserved_seqs = None

    def __repr__(self):
        return 'Phenotype: {}\nNum Genomes: {}\nOutput_dir: {}'.format(self.phenotype,
                                                                       len(self.genomes),
                                                                       self.output_dir)

    def set_phenotype(self, phenotype):
        if phenotype:
            return phenotype
        else:
            return os.path.basename(self.genome_dir)

    def add_genomes_from_dir(self, genome_dir, pheno=None):
        '''
        Takes a directory containing genome files of one phenotype and adds those
        files as genome objects to this phenotype instance's genomes list.

        param: genome_dir: String; path to the genome directory
        param: pheno: String; If not None, use this string as the phenotype 
        '''
        genomes = []
        # deciede what will be used for the phenotype name
        if pheno == None:
            self.phenotype = os.path.basename(genome_dir)
        else:
            self.phenotype = pheno

        if genome_dir:  # not none
            genomes_paths = [os.path.join(genome_dir, genome)
                             for genome in os.listdir(genome_dir)]
            for genome_path in genomes_paths:
                genomes.append(Genome(genome_path, pheno, self.output_dir))

        return genomes

    def get_conserved_sequences(self, cdhit_exec='cdhit', s='0.90', con=0.75):
        '''
        WIP
        
        TODO: Test chhit call is working as expected and add parseing of clstr
        file. Then need to find conserved seqs within the file and seperate
        those out into a new file for comparison with the other phenotype.
        
        
        Using some metrics and methods pulls out the conserved non-coding
        sequences from a collection of genomes stored in this phenotype
        instance.
        
        :param cdhit_exec: String. Path to cdhit executable default = cdhit
        :param s: String. Num 0-1 sets min length difference between rep seq \
        and subject seqs
        :param con: Int. Percentage of genomes that must have a peptide in a \
        cluster for that cluster to be considered conserved.
        '''
        
        phenotype_peptides = os.path.join(
            self.output_dir, f'{self.phenotype}_peptides.fasta')
        genome_peptides = [genome.non_coding_file for genome in self.genomes]
        # get all of the peptide file paths in one list

        cat_cmd = ['cat'] + genome_peptides + ['>', phenotype_peptides]
        cd_hit_cmd = [cdhit_exec, '-in',
                      phenotype_peptides, '-o', phenotype_peptides, '-d', '0',
                      '-p', '1', '-s', s]

        # concat all individual peptide files
        cat_call = subprocess.call(cat_cmd)
        cdhit_call = subprocess.call(cd_hit_cmd)  # run cd-hit on cated file
        
        clusters = parse_cdhit_output_file(phenotype_peptides + '.clstr')
        conserved_records = []
        
        for cluster in clusters:
            participating_genomes, rep_seq = set([]), None
            for record in cluster:
                if record[-2] == '*':
                    rep_seq = record
                participating_genomes.add(record[-1])  # add genome id
            if len(participating_genomes) / len(self.genomes) >= con:
                # percentage of genomes that participate in this cluster
                # is greater than or equal to con threshold
                conserved_records.append((rep_seq, len(cluster)))
        
        
        # at this point conserved records hold formated cdhit records of
        # representative sequences from clusters that pass the min genome
        # participation threshold
        
        #TODO Recover the actual peptide sequence of all records in the
        # conserved records list by doing some kind of query against the
        # genome identified by its genome id
                
    def pull_peptides(self, prokka_exec='prokka'):
        '''
        Wrapper around methods in the Genome class. This methos iterates
        through all genomes in a phenotype and uses prokka to make gene
        predictions, identifies non-coding regions and then translates
        those regions into six reading frames.
        '''
        for genome in self.genomes:
            genome.make_gene_predictions(path_to_exec=prokka_exec)
            genome.get_non_coding_regions()
            genome.translate_non_coding_seqs()
            genome.write_peptides_to_fasta_file()
