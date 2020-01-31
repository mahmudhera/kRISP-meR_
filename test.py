import sys
import os
import argparse

def rs2_score():
    model_file = 'V3_model_nopos.pickle'
    score = get_rs2_score('G'*30, model_file)
    print(score)
    score = get_rs2_score('G'*30, model_file)
    print(score)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("a", type=str, help="Positional argument")
    parser.add_argument("-P", "--PAM", type=str, help="Type in", nargs='+')
    args = parser.parse_args()
    print (args.PAM)
