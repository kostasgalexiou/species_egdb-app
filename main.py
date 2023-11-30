# This is a sample Python script.

# Press Mayús+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import emoji
import pandas as pd
import streamlit as st
from streamlit_extras.dataframe_explorer import dataframe_explorer
import base64
from supabase_conn import *
import tempfile
import os
from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastpCommandline, NcbiblastformatterCommandline
from collections import OrderedDict as od
import Bio 
import Bio.Blast
import Bio.Blast.Applications

# Initialize connection.
url: str = st.secrets['connections']['supabase']["SUPABASE_URL"]
key: str = st.secrets['connections']['supabase']["SUPABASE_KEY"]
supabase: Client = create_client(url, key)


def add_bg_from_local(image_file):
    with open(image_file, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read())

    st.markdown(
        f"""
        <style>
        .stApp {{
            background-image: url(data:image/{"png"};base64,{encoded_string.decode()});
            background-size: cover
        }}
        </style>
        """,
        unsafe_allow_html=True
    )


def main():
    choice = st.sidebar.radio("What do you want to do?",
                              ["Launch a new species interface :rocket:", "Select species :seedling:"], index=None)

    if choice == 'Launch a new species interface :rocket:':

        st.title('Launch a new species')
        new_species = st.text_input('Please insert latin name of new species')
        geneinfo = st.file_uploader('Upload a tab-delimited file with gene information', type='tsv')
        cdsfasta = st.file_uploader('Upload a CDS fasta file', type=['fa', 'fasta', 'fna'])
        proteinfasta = st.file_uploader('Upload a PROTEIN fasta file', type=['fa', 'fasta', 'fna'])

        if cdsfasta and proteinfasta and geneinfo:
            cds_list = getlist_from_fasta(cdsfasta)
            protein_list = getlist_from_fasta(proteinfasta)
            import_data(cds_list, new_species, 'cds')
            import_data(protein_list, new_species, 'protein')

            # d = view_species_data("protein_fastas", "Cannabis sativa")
            # st.write(d[:10])
            gene_list = pd.read_table(geneinfo, header=None).values.tolist()
            import_data(gene_list, new_species, 'gene')
    #
    #             # geneinfo_species = view_species_data("gene_info", new_species)
    #             # genedf = pd.DataFrame(geneinfo_species, columns=['geneid', 'chrom', 'start', 'end', 'transcripts', 'function'])
    #             #
    #             #
    #             # st.dataframe(cdsdf ,hide_index=True)
    #             # st.dataframe(protdf, hide_index = True)

    elif choice == "Select species :seedling:":
        all_species = ['Cannabis sativa', 'b', 'c']
        species = []
        for sp in all_species:
            if sp == 'Cucumis melo':
                species.append('melon :melon:')
            elif sp == 'Prunus persica':
                species.append('peach :peach:')
            elif sp == 'Malus domestica':
                species.append('apple :apple:')
            elif sp == 'Cannabis sativa':
                species.append('cannabis :herb:')
            else:
                species.append(sp)

        selected_species = st.sidebar.radio('select species', species, label_visibility="collapsed", index=None)

        species_menu = ['About', 'BLAST', 'Search', 'Information Extraction', 'JBrowse', 'Gene Expression', 'Tables',
                        'Downloads']

        if selected_species is not None:
            # add_bg_from_local('4656160.jpg')
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

                selected_database = st.selectbox('Select database', ['Cannabis sativa CS10 assembly - chromosomes',
                                                                     'Cannabis sativa CS10 assembly - proteins'], index=None)
                seg = 'no'
                dust = 'no'
                with st.expander('Advanced parameters'):
                    col1, col2 = st.columns(2)
                    col3, col4 = st.columns(2)

                    with col1:
                        perciden = st.slider('Percentage identity:', min_value=0.0, max_value=100.0, value=95.0, step=5.0)
                    with col2:
                        evalue = st.text_input('E value:')
                    with col3:
                        dict = od({2: 'blastp',3: 'blastx/tblastn/tblastx', 6: 'blastp/blastx/tblastn',
                                   7: 'blastn', 11: 'blastn/blastp', 28: 'blastn'})
                        dict2 = od({2: 'use with a PROTEIN database (recommended)', 3: 'use with a PROTEIN database',
                                    6: 'use with a PROTEIN database', 7: 'use with a NUCLEOTIDE database',
                                    11: 'use with a NUCLEOTIDE database', 28: 'use with a NUCLEOTIDE database (recommended)',})
                        wsize = st.selectbox('Word size', options=dict2.keys(), index=None, format_func=lambda x: str(x) + ":  " + dict2[x])
                    with col4:
                        lowcomplex = st.radio('Mask low-complexity sequences', options=['no', 'yes'], horizontal=True,
                                              index=None)

                        if lowcomplex == 'yes':
                            if 'protein' in selected_database:
                                seg = 'yes'
                            elif 'chromo' in selected_database:
                                dust = 'yes'

                # databases = supabase.storage.from_('blastDBs').download('Cannabis_sativa/chr02_1to10mbp.fa.fai').decode().split('\n')

                if selected_database and fasta_input and perciden and evalue and wsize and lowcomplex:
                    st.write(Bio.Blast.Applications.__file__, '$$$$$$$$$$$$$$$')
                    st.write(os.listdir('/home/adminuser/venv/lib'))

                    with st.spinner('Generating blast output...'):
                        fp = tempfile.NamedTemporaryFile('w+')
                        fp.write(fasta_input)
                        fp.read()
                        # fp.close()

                        with tempfile.TemporaryDirectory() as dirtemp:
                            for i in supabase.storage.from_('blastDBs').list('Cannabis_sativa'):
                                filename = i["name"]
                                with open(os.path.join(dirtemp, filename), 'wb+') as f:
                                    bucketfile = supabase.storage.from_('blastDBs').download(
                                        'Cannabis_sativa/%s' % filename)
                                    f.write(bucketfile)

                            if 'chromo' in selected_database:
                                fasta_dbname = [x for x in os.listdir(dirtemp.title().lower()) if x.endswith('.fa')][0]
                            else:
                                fasta_dbname = [x for x in os.listdir(dirtemp.title().lower()) if x.endswith('protein.faa')][0]

                            taska = ''
                            if 'chromosomes' in selected_database:
                                print('here !!!!!!!!!!!!!!!', wsize)
                                if wsize == 7:
                                    taska = 'blastn-short'
                                elif wsize == 28:
                                    taska = 'megablast'
                                elif wsize == 11:
                                    taska = 'blastn'

                                blast_cline = NcbiblastnCommandline(query=fp.name,
                                                                    db=os.path.join(dirtemp, fasta_dbname),
                                                                    evalue=float(evalue),
                                                                    perc_identity=perciden,
                                                                    word_size=wsize, task=taska, outfmt=11,
                                                                    dust=dust,
                                                                    soft_masking='false')
                            elif 'protein' in selected_database:
                                if wsize == 2:
                                    taska = 'blastp-short'
                                elif wsize == 3:
                                    taska = 'blastp'
                                elif wsize == 6:
                                    taska = 'blastp-fast'

                                blast_cline = NcbiblastpCommandline(query=fp.name,
                                                                    db=os.path.join(dirtemp, fasta_dbname),
                                                                    evalue=float(evalue),
                                                                    word_size=wsize, task=taska, outfmt=11,
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
                                    df2 = df2[df2['perc_iden'].astype(float)>=perciden]
                                    st.dataframe(df2)
                            with st.expander('Pairwise format'):
                                if '{}' in result:
                                    st.warning('No hits found!')
                                else:
                                    st.text(blast_pairwisereformat)

                    fp.close()
                    asn.close()
                    # st.text(result)
                    # df.columns = ['query', 'hit', 'perc_identity', 'length', 'mismatch', 'gap', 'qstart', 'qend',
                    #               'sstart', 'send', 'evalue', 'bitscore', 'qseq', 'sseq']

                    # AgGrid(
                    #     df,
                    #     gridOptions = GridOptionsBuilder.from_dataframe(df).build()
                    # )

            with tab3:
                pass
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


if __name__ == "__main__":
    main()
