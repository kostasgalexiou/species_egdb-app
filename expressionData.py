#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Created on Tue Dec  3 08:26:27 2013

@author: kalexiou
"""

import pandas as pd
from collections import OrderedDict as od
import streamlit as st


def genelist_dataframe(exp_data, gene_list):
    df = pd.read_table(exp_data, header=0)

    # transpose data, in order to have tissues as index, and then passing gene names as column names
    df_transp = df.iloc[:, 1:].T
    df_transp.columns = df['Name']

    if len(gene_list) == 1:
        df_gene = pd.DataFrame(df_transp[gene_list[0]])
    else:
        df_gene = df_transp[gene_list]

    grouped = df_gene.groupby(df_gene.index)

    rep2tissue = od()
    for c, gr in grouped:
        rep2tissue[c] = c.split('_rep')[0]

    df_gene['Tissue'] = [rep2tissue[x] for x in df_gene.index.values.tolist()]


    return df_gene


def genelist_dataframe2(exp_data, gene_list):
    df = pd.read_table(exp_data, header=0)
    df_transp = df.iloc[:, 1:].T
    df_transp.columns = df['Name']

    with st.expander('Average values'):
        if len(gene_list) == 1:
            df_gene2 = pd.DataFrame(df_transp[gene_list[0]])
        else:
            df_gene2 = df_transp[gene_list]

    # grouped = df_gene.groupby(df_gene.index)

    return df_gene2
