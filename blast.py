#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Created on Tue Dec  3 08:26:27 2013

@author: kalexiou
"""
import streamlit as st
import tempfile
import os
from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastpCommandline, NcbiblastformatterCommandline
from supabase import create_client, Client
import pandas as pd


# Initialize connection.
url: str = st.secrets['connections']['supabase']["SUPABASE_URL"]
key: str = st.secrets['connections']['supabase']["SUPABASE_KEY"]
supabase: Client = create_client(url, key)


def create_tmp_files(tempdir, infile):
    with open(os.path.join(tempdir, infile), 'wb+') as f:
        bucketfile = supabase.storage.from_('blastDBs').download(
            'Cannabis_sativa/%s' % infile)
        f.write(bucketfile)


def perform_blast(fasta_seq, blast_db, wsize, evalue, perciden, dust, seg):
    
    with st.spinner('Running blast...'):

        fp = tempfile.NamedTemporaryFile('w+')
        fp.write(fasta_seq)
        fp.read()

        with tempfile.TemporaryDirectory() as dirtemp:

            for i in supabase.storage.from_('blastDBs').list('Cannabis_sativa'):
                filename = i["name"]
                if 'chromo' in blast_db and 'chr' in filename:
                    file_ = filename
                    create_tmp_files(dirtemp, file_)
                elif 'prot' in blast_db and 'protein' in filename:
                    file_ = filename
                    create_tmp_files(dirtemp, file_)

            if 'chromo' in blast_db:
                fasta_dbname = [x for x in os.listdir(dirtemp.title().lower()) if x.endswith('.fa')][0]
            else:
                fasta_dbname = [x for x in os.listdir(dirtemp.title().lower()) if x.endswith('protein.faa')][0]

            task_ = ''
            if 'chromosomes' in blast_db:
                if wsize == 7:
                    task_ = 'blastn-short'
                elif wsize == 28:
                    task_ = 'megablast'
                elif wsize == 11:
                    task_ = 'blastn'

                blast_cline = NcbiblastnCommandline(query=fp.name,
                                                    db=os.path.join(dirtemp, fasta_dbname),
                                                    evalue=float(evalue),
                                                    perc_identity=perciden,
                                                    word_size=wsize, task=task_, outfmt=11,
                                                    dust=dust,
                                                    soft_masking='false')
            elif 'protein' in blast_db:
                if wsize == 2:
                    task_ = 'blastp-short'
                elif wsize == 3:
                    task_ = 'blastp'
                elif wsize == 6:
                    task_ = 'blastp-fast'

                blast_cline = NcbiblastpCommandline(query=fp.name,
                                                    db=os.path.join(dirtemp, fasta_dbname),
                                                    evalue=float(evalue),
                                                    word_size=wsize, task=task_, outfmt=11,
                                                    seg=seg)

            result, __ = blast_cline()
            asn = tempfile.NamedTemporaryFile('w+')
            asn.write(result)
            asn.read()

            cline2 = NcbiblastformatterCommandline(archive=asn.name, outfmt='6 std qseq sseq')
            blast_tabreformat, __ = cline2()
            cline3 = NcbiblastformatterCommandline(archive=asn.name, outfmt=0)
            blast_pairwisereformat, __ = cline3()

            st.success('Blast analysis is finished. You can see the results, if any, below:')
            with st.expander('Tabular format'):
                if not blast_tabreformat:
                    st.warning('No hits found!')
                else:
                    lists = [x.split('\t') for x in blast_tabreformat.rstrip('\n\r').split('\n')]
                    df = pd.DataFrame(lists, columns=['query', 'hit', 'perc_iden', 'aligned', 'mismatches', 'gaps',
                                                      'query_start', 'query_end', 'hit_start', 'hit_end', 'evalue',
                                                      'bitscore', 'query_seq', 'hit_seq'])
                    df2 = df.set_index('query')
                    df2 = df2[df2['perc_iden'].astype(float) >= perciden]
                    st.dataframe(df2)
            with st.expander('Pairwise format'):
                if '{}' in result:
                    st.warning('No hits found!')
                else:
                    st.text(blast_pairwisereformat)

        fp.close()
        asn.close()
