# This is a sample Python script.

# Press Mayús+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import pandas as pd
import base64
from supabase_conn import *
from collections import OrderedDict as od
from blast import perform_blast
import textwrap
import search

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
            # add_bg_from_local('4656160.jpg', )
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
                                                  index=None)

                        if lowcomplex == 'yes':
                            if 'protein' in selected_database:
                                seg = 'yes'
                            elif 'chromo' in selected_database:
                                dust = 'yes'

                # databases = supabase.storage.from_('blastDBs').download('Cannabis_sativa/chr02_1to10mbp.fa.fai').decode().split('\n')
                runblast = st.button('Run blast', 'runblast')

                if runblast:
                    perform_blast(fasta_seq=fasta_input, blast_db=selected_database, wsize=wsize, evalue=evalue,
                                  perciden=perciden, dust=dust, seg=seg)

            with tab3:
                entries = None
                col1, col2 = st.columns([0.4, 0.6])
                with col1:
                    keyword = st.text_input('Keyword search', placeholder='e.g. kinase, DNA polymerase, etc.')
                with col2:
                    idsearch = st.text_input('ID search', placeholder='e.g. LOC115694674, XM_030621757.1, XP_030477616.1')

                if keyword or idsearch:
                    if keyword:
                        entries = view_selected_data('gene_info', 'cannabis', 'function', keyword)
                    elif idsearch:
                        if idsearch.startswith('LOC'):
                            entries = view_selected_data('gene_info', 'cannabis', 'geneid', idsearch)
                        elif idsearch.startswith('XM'):
                            entries = view_selected_data('gene_info', 'cannabis', 'transcripts', idsearch)
                        elif idsearch.startswith('XM'):
                            entries = view_selected_data('gene_info', 'cannabis', 'proteins', idsearch)

                    if entries:
                        search.search_results(entries=entries)
                    else:
                        st.warning('No entries found...Please try with a different keyword or ID. \n\n If performing an "ID search", make sure that you abide to the gene and transcript/protein ID nomenclature')

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
