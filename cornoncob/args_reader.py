import argparse
import sys

from cornoncob import TEST_KILLERS, TEST_NICE

def get_args():
    '''
    Uses argparse to get any command line arguements that we might want
    to pass in to the program. Returns a parsed args object. Running the main
    module with --h will print out arguement descriptions.
    '''
    
    parser = argparse.ArgumentParser(description='Arguements for CoRNonCOB Piepline')
    parser.add_argument('-p1', help='Path to directory containing all genomes of the positive phenotype ')
    parser.add_argument('-p2', help='Path to directory containing all genome of the control (wild-type) phenotype')
    parser.add_argument('-o', default='.', help='Path to output directory')
    parser.add_argument('-t', default=2, help='Number of threads to use while running Prokka')
    parser.add_argument('-n', default='cornoncob_run', help='Run name')
    parser.add_argument('-k', default='prokka', help='Path to prakka executable if not in PATH variable')
    parser.add_argument('-test', default=False, help='If True, runs program in test mode')
    parser.add_argument('-s', default=0.90, help='Proportion of coverage of \
        subject sequence to representative sequence required for CD-HIT.\
        A value of 1 means sequences must be exact same length.\
        Values less than one allow subject sequences to be shorter than\
        representative sequences.')
    parser.add_argument('-c', default=0.85, help='Min proportion of genomes that\
        must participate in a cluster of similar peptides in order for that cluster\
        to be considered conserved within the phenotpe. A value of 1 will\
        mean all genomes must contribute a peptide to a cluster in order\
        for it to be considered conserved.')
     
    # add args and do some basic validations
    args = parser.parse_args()
    
    if not args.p1 or not args.p2:  # check to make sure p1 and p2 have values
        if args.test:
            args.p1 = TEST_KILLERS  # in test mode with no p1 or p2 use included data
            args.p2 = TEST_NICE
        else:
            print('Please supply two directories of different phenotypes')
            sys.exit(1)  # not in test mode so exit program
    
    if args.p1 == args.p2:  # phenotype dirs cannot be the same
        print('-p1 and -p2 must point to different loctions')
        sys.exit(3)
    
    if args.c == 0:
        print('Please set a participation value above 0')
        sys.exit(2) 
    # rest of args have default values

    return args
