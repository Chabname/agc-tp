#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
    """"
    Generate 
    Parameters:
        
    Returns:
        
    """

    path = isfile(amplicon_file)

    with gzip.open(path, "rt") as file:
        prot_dict = {}
        prot_id = ""
        for line in file:
            if line.startswith(">"):
                prot_id = line[1:].split()[0]
                prot_dict[prot_id] = ""
            else:
                prot_dict[prot_id] += line.strip()
        for i in prot_dict:
            sequence = prot_dict[i]
            if len(sequence) >= minseqlen:
                yield sequence




def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """"
    
    Parameters:
        
    Returns:
        
    """
    occurence = {}
    for sequence in read_fasta(amplicon_file, minseqlen):
        if sequence not in occurence:
            occurence[sequence] = 0
        occurence[sequence] += 1

    occu_sorted = {key: value for key, value in sorted(occurence.items(), key=lambda item: item[1], reverse = True)}

    for  seq, occu in occu_sorted.items():
        if occu >= mincount:
            yield [seq, occu]



def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """"""
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError("Sequence length ({}) is too short to be splitted in 4"
                         " chunk of size {}".format(len_seq, chunk_size))
    return [sequence[i:i+chunk_size] 
              for i in range(0, len_seq, chunk_size) 
                if i+chunk_size <= len_seq - 1]


def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]

def get_identity(alignment_list):
    """Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    id_nu = 0

    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)



def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """"
    
    Parameters:
        
    Returns:
        
    """
    seq_occu = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    kmer_dict = {}
    
 
    for index, item in enumerate(seq_occu):
        chimera = False
        chunk_list = get_chunks(item[0], chunk_size)
        id_mat = []
        for chunk in chunk_list:
            mates = search_mates(kmer_dict, chunk, kmer_size)
        if len(mates) >= 2:
            for _ in range(len(chunk_list)):
                line = []
                line.append(get_identity(nw.global_align(item[0], mates[0][0])))
                line.append(get_identity(nw.global_align(item[0], mates[1][0])))
                id_mat.append(line)
            chimera = detect_chimera(id_mat)
        kmer_dict = get_unique_kmer(kmer_dict, item[0], index, kmer_size)
        if not chimera :
            yield([item[0], item[1]])





def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """"
    
    Parameters:
        
    Returns:
        
    """
    otu = []
    seq_occu = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    for index, item  in enumerate(seq_occu):
        if index == 0:
            otu.append(item)
        else :
            for index_otu in range (len(otu)):
                align = nw.global_align(otu[index_otu][0], seq_occu[index][0], 
                        gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
                if get_identity(align) <= 97:
                    otu.append(item)
    return otu



def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    """"
    
    Parameters:
        
    Returns:
        
    """
    with open(output_file, "w", encoding="utf-8") as file:
        for index, otu in enumerate(OTU_list):
            file.write(">OTU_" + str(index +1) + " occurrence:" + str(otu[1]) + "\n")
            file.write(fill(otu[0]) + "\n")


def get_unique_kmer(kmer_dict, seq, seq_id, kmer_size):
    """"
    
    Parameters:
        
    Returns:
        
    """
    kmer_list = list(cut_kmer(seq, kmer_size))
    for kmer in kmer_list:
        if kmer in kmer_dict:
            list_id = kmer_dict[kmer]
            list_id.append(seq_id)
            kmer_dict[kmer] = list_id
        else:
            kmer_dict[kmer] = [seq_id]
    return kmer_dict



def search_mates(kmer_dict, chunk, kmer_size):
    """"
    
    Parameters:
        
    Returns:
    
    """
    kmer_list = list(cut_kmer(chunk, kmer_size))
    id_list = []
    best_match = []
    for kmer in kmer_list:
        if kmer in kmer_dict.keys():
            id_list.extend(kmer_dict[kmer])
    
    for best_id in Counter(id_list).most_common(2) :
        best_match.append(best_id[0])
    return best_match
    


def detect_chimera(perc_identity_matrix):
    """"
    
    Parameters:
        
    Returns:
    
    """
    std_list = []
    bool_result = False
    boo_seq1 = False
    bool_seq2 = False

    for perc_id in perc_identity_matrix:
        std_list.append(std(perc_id))
        if perc_id[0] > perc_id[1]:
            boo_seq1 = True
        else:
            bool_seq2 = True
        if(boo_seq1 & bool_seq2):
            break
    
    if statistics.mean(std_list) > 5:
        if(boo_seq1 & bool_seq2):
            bool_result = True
    else:
        bool_result = False
    return bool_result



def std(data):
    """
    Standard deviation

    Parameters:
        data: List of values
    Returns:
        Return the standard deviation
    """
    return statistics.stdev(data)
        

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici

#    read_fasta("data/amplicon.fasta.gz",50)


#if __name__ == '__main__':
#    main()


def main_test():
    chimerafree = chimera_removal(os.path.abspath(os.path.join("tests/test_sequences.fasta.gz")),
        200, 3, 50, 8)
    
main_test()