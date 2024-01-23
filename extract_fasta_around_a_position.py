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


def infile2dict(file_contents):

    chr2info = od()
    for line in file_contents:

        if line.startswith('#'):
            continue

        line = line.rstrip('\n\r').split('\t')

        id_, chrom, pos = line[:3]

        if chrom not in chr2info.keys():
            pos2rest = od()

        pos2rest[pos] = line[3:]
        chr2info[chrom] = pos2rest

    return chr2info


def extract_fasta(fasta, range, allele_info, output, infile_dict=od()):

    with open(op.join(os.getcwd(), output), 'w') as out:

        if allele_info:
            out.writelines('\t'.join(('#id', 'chrom', 'genome_ref_allele', 'pos', 'ref', 'alt', 'sequence\n')))
        else:
            out.writelines('\t'.join(('#id', 'chrom', 'genome_ref_allele', 'pos', 'sequence\n')))

        for k, pos_dict in infile_dict.items():

            for position, rest in pos_dict.items():

                        fasta_start_coord = int(position)-range
                        fasta_end_coord = int(position)+range
                        command = ' '.join(('samtools', 'faidx', fasta,
                                            '{0}:{1}-{2}'.format(k, fasta_start_coord, fasta_end_coord)))

                        c = subp.Popen(command, stdout=subp.PIPE, shell=True)

                        single_line = ''

                        for a in c.stdout.readlines():
                            a = a.decode("utf-8")
                            if a.startswith('>'):
                                continue
                            single_line += a.strip()
                        st.write(single_line, '$$$$$$$$$$$')
                        l_seq, genome_ref, r_seq = [single_line[:range].upper(),
                                                    single_line[range].upper(),
                                                    single_line[range+1:].upper()]

                        if allele_info == 'yes':
                            ref = rest[0]
                            alt = rest[1]

                            if len(ref) == len(alt) and len(ref) == 1:  # it is a snp
                                out.writelines('\t'.join((k + "_" + position, k, genome_ref, position, '\t'.join(rest),
                                                          l_seq + '[' + ref + '/' + alt + ']' + r_seq + '\n')))
                            else:
                                if (len(ref) != 1 or len(alt) != 1) and len(ref) > len(alt):  # it is a deletion
                                    right_seq = r_seq[len(ref[1:]):]
                                elif (len(ref) != 1 or len(alt) != 1) and len(ref) < len(alt):  # it is an insertion
                                    right_seq = r_seq
                                elif len(ref) > 1 and len(alt) > 1:
                                    print('Site {0}:{1} is a MNP'.format(k, position))

                                out.writelines('\t'.join((k + "_" + position, k, genome_ref, position, '\t'.join(rest),
                                                          l_seq + '[' + ref + '/' + alt + ']' + right_seq + '\n')))
                        if allele_info == 'no':
                            out.writelines('\t'.join((k + "_" + position, k, genome_ref, position, '\t'.join(rest),
                                                      l_seq + '@' + r_seq + '\n')))
    return

