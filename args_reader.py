import argparse
import sys


def get_args():
    '''
    Uses argparse to get any command line arguements that we might want
    to pass in to the program. Returns a parsed args object. Running the main
    module with --h will print out arguement descriptions.
    '''
    
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-p1', help='Path to directory containing all genomes of phenotype 1')
    parser.add_argument('-p2', help='Path to directory containing all genome of phenotype 2')
    parser.add_argument('-o', default='.', help='Path to output directory')
    parser.add_argument('-t', default=2, help='Number of threads')
    parser.add_argument('-n', default='corncob', help='Run name')
    parser.add_argument('-k', default='prakka', help='Path to prakka executable if not in PATH variable')
    # potentially remove
    
    
    # add args 
    args = parser.parse_args()
    if not args.p1 or not args.p2:
        print('Please supply two directories of different phenotypes')
        sys.exit(1)
    
    
    return args