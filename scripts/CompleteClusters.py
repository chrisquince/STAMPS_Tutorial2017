import sys, getopt
import os
import pandas as p
import numpy as np
import argparse
import math

from Bio import SeqIO

from collections import defaultdict
from collections import Counter

def rchop(thestring, ending):
  if thestring.endswith(ending):
    return thestring[:-len(ending)]
  return thestring


def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("scg_file", help="scg frequencies tsv")

    args = parser.parse_args()


    
    for line in open(args.scg_file):
        line = line.rstrip()
    
        tokens = line.split("\t")

        count=0

        for tok in tokens[3:]:
            if tok == "1":
                count=count+1
        if count/36. > 0.75:
            print "Cluster" + tokens[0]

    #import ipdb; ipdb.set_trace()
if __name__ == "__main__":
    main(sys.argv[1:])
