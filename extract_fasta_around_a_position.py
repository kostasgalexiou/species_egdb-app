# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 08:26:27 2013

@author: kalexiou
"""

import os.path as op
import os
from collections import OrderedDict as od
import subprocess as subp
import streamlit as st
import pysam
import pandas as pd


def infile2dict(file_contents):
    chr2info = od()
    for line in file_contents:

        if line.startswith('#'):
            continue

        line = line.rstrip('\n\r').split('\t')

        id_, chrom, pos = line[:3]

        if chrom not in chr2info.keys():
            pos2rest = od()

        pos2rest[pos + '#' + id_] = line[3:]
        chr2info[chrom] = pos2rest

    return chr2info


def extract_fasta(fasta, range_, allele_info, infile_dict=od()):

    output_list = list()
    out_df = None
    error_message = None


    for k, pos_dict in infile_dict.items():
        for p, rest in pos_dict.items():

            position, snpid = p.split('#')
            ref = pysam.FastaFile(fasta)
            string_coord = int(position) - 1

            fasta_start_coord = int(string_coord) - range_
            fasta_end_coord = int(string_coord) + range_

            try:
                l_seq = ref.fetch(k, fasta_start_coord,
                                  string_coord)  # using position because the right end of the string is open ( ref[start - end) )
                genome_ref = ref.fetch(k, string_coord, string_coord + 1)
                r_seq = ref.fetch(k, int(position), fasta_end_coord + 1)

                if allele_info == 'yes':
                    header = ['#id', 'chrom', 'pos', 'genome_ref_allele', 'ref', 'alt', 'sequence']

                    refallele = rest[0].upper()
                    altallele = rest[1].upper()

                    if len(refallele) == len(altallele) and len(refallele) == 1:  # it is a snp
                        output_list.append([snpid, k, position, genome_ref,refallele, altallele,
                                            l_seq + '[' + refallele + '/' + altallele + ']' + r_seq])
                    else:
                        if len(refallele) > 1 and len(altallele) > 1:
                            st.error('Site {0}:{1} is a MNP'.format(k, position))
                            continue
                        elif (len(refallele) > 1 or len(altallele) > 1) and len(refallele) > len(
                                altallele):  # it is a deletion
                            right_seq = ref.fetch(k, string_coord + len(refallele), fasta_end_coord + len(refallele))
                        elif (len(refallele) != 1 or len(altallele) != 1) and len(refallele) < len(
                                altallele):  # it is an insertion
                            right_seq = r_seq

                        output_list.append([snpid, k, position, genome_ref, refallele, altallele,
                                            l_seq + '[' + refallele + '/' + altallele + ']' + right_seq])
                if allele_info == 'no':
                    header = ['#id', 'chrom', 'pos', 'genome_ref_allele', 'sequence']
                    output_list.append([snpid, k, position, genome_ref,
                                        l_seq + '@' + r_seq])

                print(output_list, header)
                out_df = pd.DataFrame(data=output_list, columns=header)

            except KeyError:
                error_message = "It seems that the input file's chromosome or protein IDs are not found in the " \
                                "fasta file selected above. Make sure that you have selected the correct fasta file from the drop-down menu."

    return out_df, error_message
