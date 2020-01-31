import numpy as np
import argparse
from personalized.generate_personalized_target import *
import generate_all_candidates as candidate_generator
import pandas as pd
from calculate_priors import determine_priors_posteriors
import dna_jellyfish as jellyfish
from MLE import get_target_coverage_after_refining
from generate_adjacent_mers import generate_adjacent_mers
from get_cfd_score import get_score
from matplotlib import pyplot as plt

pam = "NGG"
grna_length = 20
candidate_length = 23
jf_count_file = "output/jf_binary_file.jf"
probability_table = []
max_k = -1
max_limit_count = 50
hist_output = 'output/k_spectrum_histo_data'
read_coverage = -1
target_coverage = -1


def parse_arguments():
    """
    parses the arguments
    :return: the parsed argument
    """
    print ('parsing the passed arguments...')
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--preprocess",
                        help="performs the initial processing on the sequencing reads and helps you find the read "
                             "coverage. If you already know the read coverage and do not wish to preprocess, "
                             "then provide the value of the read coverage with the flag -r", 
                        action="store_true")
    parser.add_argument("-e", "--em", help="performs expectation maximization to determine prior probabilities",
                        action="store_true")
    parser.add_argument("-s", "--stop", help="identify stop codons in guides and exclude those guides",
                        action="store_true")
    parser.add_argument("-v", "--detect_variant", help="uses bowtie2, samtools and pilon to detect genomic variation",
                        action="store_true")
    parser.add_argument("-t", "--test", help="test the pipeline with analysis made with reference genome",
                        action="store_true")
    parser.add_argument("-r", "--readcoverage", type=int, help="enter the read coverage as an INTEGER")
    parser.add_argument("reads_file", type=str, help="provide the filename of the reads file with path")
    parser.add_argument("target_file", type=str,
                        help="provide the filename of the file where there is the target region. Specify with full "
                             "path.") 
    parser.add_argument("-g", "--genome", type=str, help="fasta filename of genome, only needed if you want to test")
    parser.add_argument("scores_file", type=str, help="provide the name of the file where you wish to write the scores")
    parser.add_argument("max_hd", type=int,
                        help="provide the maximum hamming distance between the candidate gRNA and the genomic k-mer "
                             "that you wish for the pipeline to analyze. Make sure the hamming distance is between 0 "
                             "and 3, inclusive") 
    args_after_parsing = parser.parse_args()
    if args_after_parsing.max_hd > 3 or args_after_parsing.max_hd < 0:
        print ('Please enter correct value of hamming distance.')
        exit(-1)
    if args_after_parsing.preprocess is False and args_after_parsing.readcoverage is None:
        print ('Please either choose to preprocess or manually enter read coverage.')
        exit(-1)
    print ('Finished parsing the arguments.\n')
    return args_after_parsing


def initial_jellyfish_run(reads_file_for_jellyfish):
    """
    k-mer counting using jellyfish. takes the reads file, counts k-mers and stores in another file named
    'jf_jf_mer_count_file'
    :param reads_file_for_jellyfish: the reads file, fasta or fastq :return: None
    """
    jf_command = "jellyfish count -m " + str(
        candidate_length) + " -s 100M -o " + jf_count_file + " -t 20 -C " + reads_file_for_jellyfish
    jf_command_args = jf_command.split(" ")
    # todo: uncomment this later
    #subprocess.call(jf_command_args)
    return jf_count_file


def generate_k_spectrum_histogram(jellyfish_file, histo_output_file=hist_output):
    """
    generate the histogram using jellyfish command, then store the data in a dictionary and return that
    :param jellyfish_file: the jf binary file path
    :param histo_output_file: the file where you want to write the histo data
    :return: a dictionary containing the histogram information
    """
    histo_command = 'jellyfish histo ' + jellyfish_file
    histo_command_args = histo_command.split(' ')
    res = subprocess.check_output(histo_command_args)
    with open(histo_output_file, 'w') as f:
        f.write(res)
    histo_dataframe = pd.read_csv(histo_output_file, delimiter=' ', names=['index', 'value'])
    dic = pd.Series(histo_dataframe.value.values, index=histo_dataframe.index).to_dict()
    return dic


def plot_histogram(histo_data):
    """
    plot the k-spectrum histogram, the data is in the histo_data dictionary. Shows the image
    :param histo_data: a dictionary containing the histogram data
    :return: None
    """
    lists = sorted(histo_data.items())
    x, y = zip(*lists)
    plt.plot(x, y)
    plt.ylim(0,100)
    plt.xlim(0,20)
    plt.show()


def generate_k_spectrum_of_target_and_count(target_string, jellyfish_count_file, max_k_limit=200, k=15):
    """
    k-spectrum of target, then count the k-mers found within the target, then generate the histogram
    :type max_k_limit: int
    :param target_string: the target string
    :param k: value of k
    :param jellyfish_count_file: jellyfish binary file name
    :param max_k_limit: max value upto which the histogram is to be generated
    :return: the histogram data in a dictionary
    """
    target = target_string
    length = len(target)
    a = set()
    for i in range(length - k):
        a.add(target[i:i + k])
    lst = []
    qf = jellyfish.QueryMerFile(jellyfish_count_file)
    for substr in a:
        mer = jellyfish.MerDNA(substr)
        count = qf[mer]
        lst.append(count)
    dic = {}
    for i in range(max_k_limit):
        dic[i + 1] = lst.count(i + 1)
    return dic


def sort_second(val):
    return val[1]


def complement(seq):
    complement_char = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    bases = [complement_char[base] for base in bases]
    return ''.join(bases)


def reverse_complement(s):
    return complement(s[::-1])


def get_poisson(count, k):
    logarithm = -read_coverage * count + k * np.log(read_coverage * count)
    for i in range(2, k + 1):
        logarithm -= np.log(i)
    return np.exp(logarithm)


def get_probability(count, k):
    if k > max_k:
        return 0.0
    if probability_table[count][k] != -1:
        return probability_table[count][k]
    #total = 0.0
    for k1 in range(max_k):
        probability_table[count][k1] = get_poisson(count, k1)
        # total = total + probability_table[count][k1]
    # for k1 in range(max_k):
    #	probability_table[count][k1] = probability_table[count][k1]/total
    return probability_table[count][k]


def annotate_guides_with_score(candidates_count_dictionary, jellyfish_filename, priors, posteriors, max_hd,
                               target_string, target_coverage):
    iteration_count = 0
    list_candidates = []
    for candidate in list(candidates_count_dictionary.keys()):
        strand_type = candidates_count_dictionary[candidate]
        trie = generate_adjacent_mers(candidate, max_hd)
        value1 = value2 = 0.0
        print('processing candidate ' + candidate)
        flag = True
        for mer in trie.keys():
            if strand_type == '+':
                cp = get_score(candidate, mer)
            else:
                cp = get_score(reverse_complement(candidate), reverse_complement(mer))
            qf = jellyfish.QueryMerFile(jellyfish_filename)
            merDNA = jellyfish.MerDNA(mer)
            rev_comp_merDNA = jellyfish.MerDNA(reverse_complement(mer))
            k = max(qf[merDNA], qf[rev_comp_merDNA])
            if k <= 0:
                continue
            if k >= max_k:
                flag = False
                break
            p = float(target_string.count(mer))
            accum = 0.0
            for count in range(1, max_limit_count):
                probability = get_probability(count, k)
                p_count = priors[count]
                p_k = posteriors[k]
                new_val = 1.0 * probability * count * p_count / p_k
                accum = accum + new_val
            value1 = value1 + cp * p
            value2 = value2 + cp * accum
        if value1 <= 0.0 or flag is False:
            continue
        score = 1.0 * value2 / (value1 * target_coverage)
        qf = jellyfish.QueryMerFile(jellyfish_filename)
        merDNA = jellyfish.MerDNA(candidate)
        k = max(qf[merDNA], qf[jellyfish.MerDNA(reverse_complement(candidate))])
        list_candidates.append((candidate, score, k, trie, strand_type))
        iteration_count = iteration_count + 1
        print('processed ' + str(iteration_count) + 'th gRNA: ' + candidate + ' with score= ' + str(score))
    print('DONE! Sorting...')
    list_candidates.sort(key=sort_second)
    print('Final list:')
    f = open('scores', 'w')
    for annotated_candidate in list_candidates:
        print(annotated_candidate)
        f.write(str(annotated_candidate[1]) + '\n')
    f.close()
    return list_candidates


def krispmer_main(parsed_args):
    """
    takes the arguments parsed by the arg parser, then performs the pipeline of the tool
    :param parsed_args: the arguments passed to the programs
    :return: a list of all guides annotated with off-target activity
    """
    # get all arguments
    # already have in this version

    # initial jellyfish run
    print('doing an initial run of Jellyfish to count k-mers in reads. Please wait...')
    jellyfish_binary_file = initial_jellyfish_run(parsed_args.reads_file)
    print('Completed the initial run.\n')

    # personalized gRNA as a string
    if args.detect_variant:
        print ('generating personalized version of the target\n')
        modified_target_string = detect_variant(parsed_args.target_file, parsed_args.reads_file)
        print ('personalized target identified\n')
    else:
        modified_target_string = read_target_region(parsed_args.target_file)

    # generate k-mer spectrum histogram data
    print('generating histogram data from the initial Jellyfish database.')
    histogram_data_dictionary = generate_k_spectrum_histogram(jellyfish_binary_file)
    print('finished generating histogram data\n')

    # do the preprocess
    global read_coverage
    if parsed_args.preprocess:
        print ('running the initial processing.')
        plot_histogram(histogram_data_dictionary)
        print ('input the read coverage...')
        read_coverage = float(input())  # type: float
    else:
        read_coverage = parsed_args.readcoverage

    # determine all candidate list
    print('generating list of potential candidates...')
    candidates_count_dictionary = candidate_generator.get_list_of_candidates(modified_target_string, pam, grna_length, args.stop)
    print('finished generating list of potential candidates...\n')

    # determine priors, posteriors and read-coverage using EM
    global max_k
    print('calculating priors...')
    priors, posteriors, read_coverage, inversion_point = determine_priors_posteriors(histogram_data_dictionary,
                                                                                     parsed_args.em,
                                                                                     max_limit_count)
    max_k = len(posteriors)
    print('finished calculating priors and posteriors\n')

    # initialize poisson probability table
    global probability_table
    probability_table = [[-1] * max_k for i in range(max_limit_count)]  # type: List[List[int]]

    # perform MLE to determine target coverage
    global target_coverage
    k_spectrum_data_in_target = generate_k_spectrum_of_target_and_count(modified_target_string, jellyfish_binary_file)
    target_coverage = get_target_coverage_after_refining(k_spectrum_data_in_target, read_coverage, inversion_point)

    # annotate all guides
    print('processing total ' + str(len(list(candidates_count_dictionary.keys()))) + ' candidate gRNAs')
    list_candidates = annotate_guides_with_score(candidates_count_dictionary, jellyfish_binary_file, priors, posteriors,
                                                 parsed_args.max_hd, modified_target_string, target_coverage)
    return list_candidates


if __name__ == '__main__':
    args = parse_arguments()
    gRNAs = krispmer_main(args)
    print (gRNAs)
