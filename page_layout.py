#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Created on Tue Dec  3 08:26:27 2013

@author: kalexiou
"""

import streamlit as st
from blast import perform_blast
from collections import OrderedDict as od
from supabase_conn import *
import search


if 'keyword' not in st.session_state or 'idsearch' not in st.session_state or \
    'runblast' not in st.session_state:
    st.session_state.keyword = ''
    st.session_state.idsearch = ''
    st.session_state.runblast = ''


def submit():
    st.session_state.keyword = st.session_state.keywidget
    st.session_state.idsearch = st.session_state.idwidget
    st.session_state.runblast = st.session_state.blastwidget
    st.session_state.keywidget = ''
    st.session_state.idwidget = ''
    st.session_state.fastawidget = ''


def generate_page(species):

    species_menu = ['About', 'BLAST', 'Search', 'Information Extraction', 'JBrowse', 'Gene Expression', 'Tables',
                    'Downloads']

    tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8 = st.tabs(species_menu)

    with tab1:
        pass
    with tab2:
        fasta_input = st.text_area('Please paste fasta sequence:', height=300,
                                   placeholder='>input_fasta\nATGGTTCATTGTGATGATGATAATCCTCAT'
                                               'GGTGTGGATGGTTTGGTTGCCAAAATTGCTGATTTTGGTCTCTC'
                                               'TAAGAGAAGTAATGATGATTTTGAGAAGCGAGTGAGAGGTACTT'
                                               'CGTTGTATCTGTCTCCGAGAGTGTTTCTTGAAGAAACTCAAGAC'
                                               'CAGGCCTCTGACATTTGGGCTCTTGGGTGTGTTGTTCTTGAGAT'
                                               'GTTGACTGGAAGTACGCCATGGGATGTAAATTTGGGGCAGGAAG'
                                               'AGCTGTTCAGCAAGATTTCTACTGAAACACCTTCTCTACCATCT')

        # fasta_input = st.session_state.fasta_input

        selected_database = st.selectbox('Select database', ['Cannabis sativa CS10 assembly - chromosomes',
                                                             'Cannabis sativa CS10 assembly - proteins'],
                                         index=None)
        seg = 'no'
        dust = 'no'
        with st.expander('Advanced parameters'):
            col1, col2 = st.columns(2)
            col3, col4 = st.columns(2)

            with col1:
                perciden = st.slider('Percentage identity:', min_value=0.0, max_value=100.0, value=95.0,
                                     step=5.0)
            with col2:
                evalue = st.text_input('E value:', placeholder='0.0001, 1e-100', value='1e-50')
            with col3:
                dict = od({2: 'blastp', 3: 'blastx/tblastn/tblastx', 6: 'blastp/blastx/tblastn',
                           7: 'blastn', 11: 'blastn/blastp', 28: 'blastn'})
                dict2 = od({2: 'use with a PROTEIN database (default)', 3: 'use with a PROTEIN database',
                            6: 'use with a PROTEIN database', 7: 'use with a NUCLEOTIDE database',
                            11: 'use with a NUCLEOTIDE database', 28: 'use with a NUCLEOTIDE database '
                                                                      '(default)', })

                if not selected_database:
                    wsize = st.selectbox('Word size', options=dict2.keys(), index=None,
                                         format_func=lambda x: str(x) + ":  " + dict2[x])
                if selected_database and 'chrom' in selected_database:
                    wsize = st.selectbox('Word size', options=dict2.keys(), index=5,
                                         format_func=lambda x: str(x) + ":  " + dict2[x])
                elif selected_database and 'protein' in selected_database:
                    wsize = st.selectbox('Word size', options=dict2.keys(), index=0,
                                         format_func=lambda x: str(x) + ":  " + dict2[x])
                with col4:
                    lowcomplex = st.radio('Mask low-complexity sequences', options=['no', 'yes'],
                                          horizontal=True,
                                          index=0)

                if lowcomplex == 'yes':
                    if 'protein' in selected_database:
                        seg = 'yes'
                    elif 'chromo' in selected_database:
                        dust = 'yes'

        # databases = supabase.storage.from_('blastDBs').download('Cannabis_sativa/chr02_1to10mbp.fa.fai').decode().split('\n')
        st.button('Run blast', key='blastwidget', on_click=submit)
        runblast = st.session_state.runblast

        if runblast:
            perform_blast(fasta_seq=fasta_input, blast_db=selected_database, wsize=wsize, evalue=evalue,
                          perciden=perciden, dust=dust, seg=seg)

        if st.session_state['runblast'] != '':
            st.session_state['runblast'] = ''

    with tab3:
        entries = None
        col1, col2 = st.columns([0.4, 0.6])
        with col1:
            st.text_input('Keyword search', placeholder='e.g. kinase, DNA polymerase, etc.', on_change=submit,
                          key='keywidget')
            keyword = st.session_state.keyword
        with col2:
            st.text_input('ID search', placeholder='e.g. LOC115694674, XM_030621757.1, XP_030477616.1',
                          on_change=submit, key='idwidget')
            idsearch = st.session_state.idsearch

        if keyword and not idsearch:
            entries = view_selected_data('gene_info', 'cannabis', 'function', keyword)

            if entries:
                search.search_results(entries=entries, term=keyword)
                keyword = ''
            else:
                st.warning(
                    'No entries found...Please try with a different keyword')

        if idsearch and not keyword:
            if idsearch.startswith('LOC'):
                entries = view_selected_data('gene_info', 'cannabis', 'geneid', idsearch)
            elif idsearch.startswith('XM'):
                entries = view_selected_data('gene_info', 'cannabis', 'transcripts', idsearch)
            elif idsearch.startswith('XM'):
                entries = view_selected_data('gene_info', 'cannabis', 'proteins', idsearch)

            if entries:
                search.search_results(entries=entries, term=idsearch)
            else:
                st.warning(
                    'No entries found...Please try with a different ID. \n\n '
                    'Make sure that you abide to the gene and transcript/protein ID nomenclature')

        elif keyword and idsearch:
            st.warning('Please use only one of the search fields at a time.')

        if st.session_state['keyword'] != '' or st.session_state.idsearch != '':
            st.session_state['keyword'] = ''
            st.session_state.idsearch = ''

    with tab4:
        pass
    with tab5:
        pass
    with tab6:
        pass
    with tab7:
        pass
    with tab8:
        pass