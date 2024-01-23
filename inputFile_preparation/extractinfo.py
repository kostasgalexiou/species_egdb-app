#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Created on Tue Dec  3 08:26:27 2013

@author: kalexiou
"""

import gzip
from collections import OrderedDict as od

gene2info = od()

with gzip.open('GCF_900626175.2_cs10_genomic_no-comment-lines.gff.gz', 'rt') as f:
#with gzip.open('test.gz', 'rt') as f:
    f1 = f.readlines()
    pseudogene = list()
    geneinfo = list()
    gname = ''
    gname_set = set()
    for n in range(len(f1)):
        if f1[n].startswith('#'):
            continue
        line = f1[n].rstrip('\n').split('\t')
        annotation, rest, chrom, start, end = line[2], line[8], line[0], line[3], line[4]

        if chrom == 'Cs10.MT' or annotation == 'region':
            continue

        if n < len(f1) - 1:
            next_line = f1[n + 1].rstrip('\n').split('\t')
            next_annotation = next_line[2]

        if next_annotation == 'gene' or n == len(f1) -1:
            gene2info[gname] = geneinfo
            geneinfo = list()

        if annotation in ['CDS', 'exon', 'cDNA_match']:
            continue

        if annotation == 'pseudogene':
            pseudog = rest.split(';',1)[0].split('-')[1]
            pseudogene.append(pseudog)

        elif annotation=="gene":
            geneinfo.append(chrom)
            geneinfo.append(start)
            geneinfo.append(end)

            for x in rest.split(';'):
                if x.startswith('Name'):
                    gname = x.split('=')[1]

        elif rest.split(';')[1].split('-')[1] in pseudogene:
                continue
        else:
            rnaid = rest.split(';',1)[0].split('-')[1]

            if 'product' in rest:
                index = [n for n,x in enumerate(rest.split(';')) if x.startswith('product')][0]
                function = rest.split(';')[index].split('=')[1]
            else:
                function = '-'
            geneinfo.append([rnaid, function])


with open('gene_info.tab', 'w') as out:
    for k, v in gene2info.items():
        if k in ['rps4', 'sdh3']:
            continue
        print(k, v)
        chromo, start, end = v[:3]
        allele_lists = v[3:]
        allele_set = set([''.join(x[1].split(' ')[:-1]) for x in allele_lists])

        if len(allele_lists) > 1:
            alleles = ';'.join([x[0] for x in allele_lists])
            if len(allele_set) == 1:
                function = allele_lists[0][1].rsplit(' ',1)[0]
            else:
                function = ';'.join([x[1] for x in allele_lists])
        else:
            alleles, function = allele_lists[0]

        out.writelines('\t'.join((k, chromo, str(start), str(end), alleles, function+'\n')))

#print (gene2info['LOC115713840'])