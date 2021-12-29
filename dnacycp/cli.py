"""DNAcycP - DNA Sequence Cyclizability Prediction

Usage:
    dnacycp-cli fasta <inputfile> <basename>
    dnacycp-cli txt <inputfile> <basename>
    dnacycp-cli (-h | --help)
    
Arguments:
    <inputfile> Input file name.
    <basename>  Output file name base.
    
Options:
    -h --help   Show this screen.
    
"""
from docopt import docopt
from dnacycp import cycle_fasta, cycle_txt
import keras
import pandas as pd
import numpy as np
from numpy import array
from Bio import SeqIO

def main():
    arguments = docopt(__doc__)
    print("Input file: "+arguments['<inputfile>'])

    if arguments['fasta']:
        cycle_fasta(arguments['<inputfile>'],
            arguments['<basename>'])
    elif arguments['txt']:
        cycle_txt(arguments['<inputfile>'],
            arguments['<basename>'])
            
if __name__ == "__main__":
    main()
