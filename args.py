import argparse
import sys


def get_args():
    '''
    Uses argparse to get any command line arguements that we might want
    to pass in to the program. Returns a parsed args object
    '''
    parser = argparse.ArgumentError()
    parser.add_arguement('p1', help='Path to directory containing all genomes of phenotype 1')
    parser.add_arguement('p2', help='Path to directory containing all genome of phenotype 2')
    parser.add_arguement('-o', default='.', help='Path to output directory')
    parser.add_arguement('-n', default='corncob', help='Run name')
    parser.add_arguement('-k', default='prakka', help='Path to prakka executable if not in PATH variable')
    
    
    # add args 
    parser.parse_args()
    if not parser.p1 or not parser.p2:
        print('Please supply two directories of different phenotypes')
        sys.exit(1)
    
    
    return parser