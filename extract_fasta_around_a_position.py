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


def infile2dict(file_contents):
    chr2info = od()
    for line in file_contents:

        if line.startswith('#'):
            continue

        line = line.rstrip('\n\r').split('\t')

        id_, chrom, pos = line[:3]

        if chrom not in chr2info.keys():
            pos2rest = od()

        pos2rest[pos+'#'+id_] = line[3:]
        chr2info[chrom] = pos2rest

    return chr2info


def extract_fasta(fasta, range, allele_info, output, infile_dict=od()):
    with open(op.join(os.getcwd(), output), 'w') as out:

        if allele_info:
            out.writelines('\t'.join(('#id', 'chrom', 'genome_ref_allele', 'pos', 'ref', 'alt', 'sequence\n')))
        else:
            out.writelines('\t'.join(('#id', 'chrom', 'genome_ref_allele', 'pos', 'sequence\n')))

        for k, pos_dict in infile_dict.items():
            for p, rest in pos_dict.items():

                position, snpid = p.split('#')
                ref = pysam.FastaFile(fasta)
                string_coord = int(position) - 1

                fasta_start_coord = int(string_coord) - range
                fasta_end_coord = int(string_coord) + range

                l_seq = ref.fetch(k, fasta_start_coord, string_coord)  # using position because the right end of the string is open ( ref[start - end) )
                genome_ref = ref.fetch(k, string_coord, string_coord+1)
                r_seq = ref.fetch(k, int(position), fasta_end_coord+1)

                if allele_info == 'yes':
                    refallele = rest[0].upper()
                    altallele = rest[1].upper()

                    if len(refallele) == len(altallele) and len(refallele) == 1:  # it is a snp
                        out.writelines('\t'.join((snpid, k, position, genome_ref, '\t'.join(rest),
                                                  l_seq + '[' + refallele + '/' + altallele + ']' + r_seq + '\n')))
                    else:
                        if len(refallele) > 1 and len(altallele) > 1:
                            st.error('Site {0}:{1} is a MNP'.format(k, position))
                            continue
                        elif (len(refallele) > 1 or len(altallele) > 1) and len(refallele) > len(altallele):  # it is a deletion
                            right_seq = ref.fetch(k, string_coord+len(refallele), fasta_end_coord+len(refallele))
                        elif (len(refallele) != 1 or len(altallele) != 1) and len(refallele) < len(altallele):  # it is an insertion
                            right_seq = r_seq

                        out.writelines('\t'.join((snpid, k, position, genome_ref, '\t'.join(rest),
                                                  l_seq + '[' + refallele + '/' + altallele + ']' + right_seq + '\n')))
                if allele_info == 'no':
                    out.writelines('\t'.join((snpid, k, position, genome_ref, '\t'.join(rest),
                                              l_seq + '@' + r_seq + '\n')))

    return
