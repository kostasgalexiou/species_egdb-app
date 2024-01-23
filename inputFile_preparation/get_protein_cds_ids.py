#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Created on Tue Dec  3 08:26:27 2013

@author: kalexiou
"""

import gzip

geneset = set()
with open('gene_info.tsv') as f:
    for line in f:
        gid = line.split('\t')[0]
        geneset.add(gid)

with open('protinfo.tab', 'w') as out:
    with gzip.open('GCF_900626175.2_cs10_genomic_no-comment-lines.gff.gz', 'tr') as f:
        for line in f:
            line = line.rstrip('\n\r').split('\t')

            if line[2] != 'CDS':
                continue

            geneid = [x.split('=')[1] for x in line[-1].split(';') if x.startswith('gene=')][0]

            if geneid not in geneset:
                continue

            proteinid = [x.split('=')[1] for x in line[-1].split(';') if x.startswith('Name=XP')][0]
           # print(geneid, proteinid)

from collections import OrderedDict as od
gid2prot = od()
with open('df2', 'w') as out:
    with open('df') as f:
        for line in f:
            a, b = line.rstrip('\n\r').split('\t')
            if a not in gid2prot.keys():
                protlist = []

            protlist.append(b)
            gid2prot[a] = protlist

    for g in list(geneset):
        if g in gid2prot.keys():
            out.writelines('\t'.join((g, ';'.join(gid2prot[g])+'\n')))
        else:
            out.writelines('\t'.join((g, '-\n')))
