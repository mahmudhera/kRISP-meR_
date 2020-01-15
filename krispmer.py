import numpy as np
import argparse

def parse_arguments():
    """
    parses the arguments
    :return: the parsed argument
    """
    global reads_file, target_region_file, max_hd
    print ('parsing the passed arguments...')
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--preprocess",
                        help="performs the initial processing on the sequencing reads and helps you find the read coverage. If you already know the read coverage and do not wish to preprocess, then provide the value of the read coverage with the flag -r",
                        action="store_true")
    parser.add_argument("-e", "--em", help="performs expectation maximization to determine prior probabilities",
                        action="store_true")
    parser.add_argument("-t", "--test", help="test the pipeline with genomic analysis",
                        action="store_true")
    parser.add_argument("reads_file", type=str, help="provide the filename of the reads file with path")
    parser.add_argument("target_file", type=str,
                        help="provide the filename of the file where there is the target region. Specify with full path.")
    parser.add_argument("-g", "--genome", type=str, help="fasta fiename of genome, only needed if you want to test")
    parser.add_argument("scores_file", type=str, help="provide the name of the file where you wish to write the scores")
    parser.add_argument("max_hd", type=int,
                        help="provide the maximum hamming distance between the candidate gRNA and the genomic k-mer that you wish for the pipeline to analyze. Make sure the hamming distance is between 0 and 3, inclusive")
    args = parser.parse_args()
    if args.max_hd > 3 or args.max_hd < 0:
        print ('Please enter correct value of hamming distance.')
        exit(-1)
    print ('Finished parsing the arguments.\n')
    return args


if __name__ == '__main__':
    print ('Should test the parser now.')
    print ('Extensively...')
