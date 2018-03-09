'''anchors_generator.py

This python file runs to read in fasta file of V genes or J genes in T-cell
or B-cell to translate DNA genes to amino acids and then finding the index of
the last occurrence of cysteine in the V genes or the first occurence of
phenylalanine followed by glycine in J genes.
'''
import argparse
import sys
from Bio import SeqIO
import csv
import math
import re

results = []
indexs = []
error_results =[]
error_indexs = []
sequence = []

def main(args):

    # Read in command line arguments to variables
    v_or_f = args.t
    input_dir = args.i
    output_dir = args.o

    if (v_or_f == "V" or v_or_f == "v"):
        parse_v_genes(input_dir)
        generate_error_file(output_dir)
    else:
        parse_j_genes(input_dir)

    generate_anchor_file(output_dir)


def parse_j_genes(infile):
    '''Find the anchors in a j genes file

    Attributes
    ----------
        infile (str): full path to input file
    '''
    for seq_record in SeqIO.parse(infile, "fasta"):
        results.append(seq_record.description)
        ind = []
        # try out 3 different frames
        for i in range(3):
            # split record to match triplets
            seq_record_temp = seq_record.seq[i:]
            floor = math.floor(len(seq_record_temp)/3)
            index = len(seq_record_temp) - (floor*3)
            while index != 0:
                seq_record_temp = seq_record_temp[:-(index)]
            # translating from dna to amino acid and find first FG
            # translating from dna to amino acid and find first (F/W)X(S/T)
            traslated_seq = seq_record_temp.translate()
            m = re.search('[FW]G.?G[ST]', traslated_seq)
            if m:
                position = m.start()
            index_F = ((position*3)+i)
            ind.append(index_F)
        # look for only positive indexs
        pos_idx = [i for i in ind if i >=0]
        if len(pos_idx) != 0:
                indexs.append(str(min(pos_idx)))
        else:
                indexs.append(str(0))


def parse_v_genes(infile):
    '''Find the anchors in a v genes file

    Attributes
    ----------
        infile (str): full path to input file
    '''
    for seq_record in SeqIO.parse("genomicVs.fasta", "fasta"):
        floor = math.floor(len(seq_record.seq)/3)
        index = len(seq_record.seq) - (floor*3)
        seq_record.seq = seq_record.seq[:-(index)]
        anchor_index = str(seq_record.seq.translate().rfind('C')*3)

        #filter out abnormal V genes
        threshold = len(seq_record.seq)/2
        if int(anchor_index) > threshold:
            results.append(seq_record.description)
            indexs.append(anchor_index)
        else:
            error_indexs.append(anchor_index)
            error_results.append(seq_record.description)
            sequence.append(seq_record.seq)


def generate_anchor_file(fileName):
    '''generate the anchor file for both v genes and j genes

    Attributes
    ----------
        fileName (str): file name for the anchor file generated
    '''
    with open(fileName, 'w') as csv_file:
        fieldnames = ['gene','anchor_index']
        csv_writer = csv.DictWriter(csv_file,fieldnames=fieldnames,delimiter=';')
        csv_writer.writeheader()

        for result, index in zip(results, indexs):
            csv_writer.writerow({'gene': result, 'anchor_index': index})


def generate_error_file(fileName):
    # only works for csv files
    '''Generate the error anchor file for v genes with no C or C apearance
       in the begenning of the chain

    Attributes
    ----------
        fileName (str): file name for the complementary anchor file generated
    '''
    fileName = fileName.split('.csv')[0] + '_error.csv'
    with open(fileName, 'w') as csv_file:
        fieldnames = ['gene','sequence','anchor_index']
        csv_writer = csv.DictWriter(csv_file,fieldnames=fieldnames,delimiter=';')
        csv_writer.writeheader()

        for error_result, seq, error_index in zip(error_results, sequence, error_indexs):
            csv_writer.writerow({'gene': error_result, 'sequence':seq,'anchor_index': error_index})

if __name__ == '__main__':

    # Set commend line arugments
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', help = 'type V or J')
    parser.add_argument('-i', help = 'path to the input file')
    parser.add_argument('-o', help = 'path to the output file')
    args = parser.parse_args()

    # message for wrong command line
    if (args.t == None or args.i == None or args.o == None):
        print("Command line arugment error\nCorrect Usage:\npython anchors_generator.py -t {V,J} -i <full path of input file> -o <full path to output file>")
        sys.exit()
    main(args)
