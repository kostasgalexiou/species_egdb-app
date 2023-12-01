#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Created on Tue Dec  3 08:26:27 2013

@author: kalexiou
"""

import streamlit as st
import pandas as pd
from awesome_table import AwesomeTable
import time
from supabase_conn import view_selected_data
import textwrap


def search_results(entries):

    time_string = time.strftime("%m/%d/%Y-%H:%M:%S", time.localtime())
    wrapper = textwrap.TextWrapper(width=70)

    st.write('\n\n')
    st.markdown("<h3 style='text-align: center; color: Green;'>Output results</h3>", unsafe_allow_html=True)
    st.write('\n\n\n')
    st.markdown("<h5 style='text-align: left; color: green;'>Gene information:</h5>", unsafe_allow_html=True)
    geneinfo = pd.DataFrame(entries)
    geneinfo.columns = ['geneid', 'chrom', 'start', 'end', 'transcript', 'protein', 'function', 'species']
    AwesomeTable(geneinfo)
    csv = geneinfo.to_csv(header=True, index=False)
    st.download_button(':blue[Download gene information]', csv,
                       file_name='gene-info_{}.csv'.format(time_string), mime='txt/csv')

    cds_ids = [x[4] for x in entries]
    protein_ids = [x[5] for x in entries]

    cds_entries = []
    for cds_id in cds_ids:
        c_list = cds_id.split(';')
        c_entries = [view_selected_data('cds_fastas', 'cannabis', 'seqid', c) for c in c_list]

        for c in c_entries:
            if c:
                cds_entries.append(['>' + c[0][0], c[0][1]])  # exclude species info

    protein_entries = []
    for p_id in protein_ids:
        p_list = p_id.split(';')
        pr_entries = [view_selected_data('protein_fastas', 'cannabis', 'seqid', p)[:2] for p in p_list]
        for p in pr_entries:
            if p:
                protein_entries.append(['>' + p[0][0], p[0][1]])  # exclude species info

    st.write('\n\n\n')

    st.markdown("<h5 style='text-align: left; color: green;'> mRNA sequence(s):</h5>", unsafe_allow_html=True)

    # string wrapping
    c_wrapped_string = ''
    for i in cds_entries:
        c_wrapped_string += i[0] + '\n'
        string = wrapper.fill(
            text=i[1])
        c_wrapped_string += string + '\n\n'
    st.code(c_wrapped_string)

    st.markdown("<h5 style='text-align: left; color: green;'> protein sequence(s):</h5>", unsafe_allow_html=True)

    # string wrapping
    p_wrapped_string = ''
    for i in protein_entries:
        p_wrapped_string += i[0] + '\n'
        string = wrapper.fill(
            text=i[1])
        p_wrapped_string += string + '\n\n'
    st.code(p_wrapped_string)

    whole_string = c_wrapped_string + '\n' + p_wrapped_string
    st.download_button(':blue[Download sequence(s)]', whole_string,
                       file_name='Fasta_sequences_{}.fa'.format(time_string))
