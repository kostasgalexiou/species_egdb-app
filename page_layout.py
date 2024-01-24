#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Created on Tue Dec  3 08:26:27 2013

@author: kalexiou
"""
import tempfile
from blast import perform_blast
from supabase_conn import *
import search
from extract_fasta_around_a_position import *
import os.path as op
import streamlit.components.v1 as components
from expressionData import *
import plotly.express as px
from collections import Counter
import numpy as np
import os

pd.options.plotting.backend = 'plotly'


def tab_config():
    st.markdown("""
                <style>
            
                    .stTabs [data-baseweb="tab-list"] {
                        gap: 2px;
                    }
            
                    .stTabs [data-baseweb="tab"] {
                        height: 30px;
                        white-space: pre-wrap;
                        background-color: #FF8080;
                        border-radius: 4px 4px 0px 0px;
                        gap: 1px;
                        padding-top: 10px;
                        padding-bottom: 10px;
                    }
            
                    .stTabs [aria-selected="true"] {
                        background-color: #FFFFFF;
                    }
            
                </style>""", unsafe_allow_html=True)


def get_avg_value_per_tissue(melted_df):  # Tissue, Replicates, Gene, Counts

    gene_grouped = melted_df.groupby('Gene')

    dict_for_df = od()

    gene_list = list()
    tissue_list = list()
    count_list = list()

    for g, gr in gene_grouped:
        tissue_counts = Counter(gr['Tissue'])

        # generate a gene list with a length equal at the number of tissues
        gene_list.extend([g] * len(tissue_counts.keys()))

        tissue_list.extend(sorted(tissue_counts.keys()))

        t2values = od()
        values = list()
        tissue_grouped = gr.groupby('Tissue')
        for t, gr2 in tissue_grouped:
            lista = gr2.Counts.values.tolist()
            values.append([min(lista), max(lista)])
            t2values[t] = values
            values = list()

            average = np.average(gr2.Counts.values.tolist())
            count_list.append(str(average))

    # generate the dictionary that will be passed into a dataframe
    dict_for_df['Gene'] = gene_list
    dict_for_df['Tissue'] = tissue_list
    dict_for_df['Average'] = count_list

    new_df = pd.DataFrame.from_dict(dict_for_df)

    return new_df


def extract_function(v_file, dirtemp, ffile, srange, ainfo):

    with open(op.join(dirtemp, ffile), 'wb+') as f:

        contents = v_file.getvalue().decode("utf-8").split('\n')[:-1]

        inf_dict = infile2dict(contents)

        bucketfile = supabase.storage.from_('blastDBs').download('Cannabis_sativa/%s' % ffile)
        f.write(bucketfile)

        dataf, error_msg = extract_fasta(fasta=op.join(dirtemp, ffile), range_=srange,
                                         allele_info=ainfo,
                                         infile_dict=inf_dict)

        st.write(dataf)

        if error_msg:
            st.warning(error_msg)
        else:
            dataf_csv = dataf.to_csv(header=True, index=False, sep='\t')

            fname = '%s_plusFasta.tab' % v_file.name.split('.')[0]

            outfile = op.join(os.getcwd(), fname)
            st.write(outfile)
            with open(outfile, 'w') as out:
                out.write(dataf_csv)


    return


def generate_page(species):
    st.header('_Cannabis sativa_ portal', divider='rainbow')
    species_menu = ['About ❓', 'BLAST 💥', 'Search 🔎', 'Information Extraction ℹ️', 'JBrowse 🧬', 'Gene Expression 📊',
                    'Tables',
                    'Downloads ⬇️']

    # tab_config()

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
                                               'AGCTGTTCAGCAAGATTTCTACTGAAACACCTTCTCTACCATCT',
                                   key='blastin_%s' % species)

        selected_database = st.selectbox('Select database', ['Cannabis sativa CS10 assembly - chromosomes',
                                                             'Cannabis sativa CS10 assembly - proteins'], index=None,
                                         key='blastdb_%s' % species)
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

        runblast = st.button('Run blast')

        if runblast:
            perform_blast(fasta_seq=fasta_input, blast_db=selected_database, wsize=wsize, evalue=evalue,
                          perciden=perciden, dust=dust, seg=seg)

    with tab3:
        entries = None
        col1, col2 = st.columns([0.4, 0.6])
        with col1:
            keyword = st.text_input('Keyword search', placeholder='e.g. kinase, DNA polymerase, etc.',
                                    key='keywidget_%s' % species.split(' ')[0])

        with col2:
            idsearch = st.text_input('ID search', placeholder='e.g. LOC115694674, XM_030621757.1, XP_030477616.1',
                                     key='idwidget_%s' % species.split(' ')[0])

        if keyword and not idsearch:
            entries = view_selected_data('gene_info', 'cannabis', 'function', keyword)

            if entries:
                search.search_results(entries=entries, term=keyword)
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
            st.warning('Please use only one search field at a time.')

    with tab4:
        variant_file = st.file_uploader('Upload SNP or INDEL file', type=['tsv', 'tab'],
                                        key='varfilewidget_%s' % species.split(' ')[0])

        if variant_file:
            fasta_names = [i["name"] for i in supabase.storage.from_('blastDBs').list('Cannabis_sativa') if
                           i["name"].endswith('.fa') or i["name"].endswith('.faa') or i["name"].endswith('.fasta')]
            fasta_file = st.selectbox('Select fasta file', options=fasta_names)

            seq_range = st.number_input(
                'Enter size of the fasta to be extracted around the variant (Default: 50nt)',
                50)
            allele_info = st.radio('Do you have REF and ALT allele information in the input file?',
                                   options=['yes', 'no'], index=0, horizontal=True)

            run_analysis = st.button('Get fasta sequences')

            if run_analysis:
                with st.spinner('Please wait...'):
                    extract_function(v_file=variant_file, dirtemp=os.getcwd(), ffile=fasta_file, srange=seq_range,
                                     ainfo=allele_info)

                # if dataf_csv:
                #     downl = st.download_button('download', data=dataf_csv, file_name=fname)
                #     if downl:
                #         st.success('Results are saved in %s' % outfile)

    with tab5:
        iframe_src = "http://localhost/jbrowse-1.16.11/?data=data%2Fjson%2F{0}%2Fcs10&loc=Cs10.Chr07%3A1..71238074&tracks=cs10%2Ccannabis-cs10_genes&highlight=".format(
            species.split(' ')[0])
        components.iframe(iframe_src, width=1000, height=500, scrolling=True)
        st.write('[Open in Jbrowse](%s)' % iframe_src)
    with tab6:
        st.markdown("<h3 style='text-align: center; color: white;'> Expression Data </h3>", unsafe_allow_html=True)
        st.markdown("<br></br>", unsafe_allow_html=True)

        col1, col2, col3 = st.columns([0.4, 0.2, 0.4])

        with col1:
            counts = st.file_uploader('$\\bold{Upload\ gene\ counts:}$', type=['csv', 'tab'],
                                      key='countwidget_%s' % species)

        with col2:
            st.markdown('$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \\bold{OR}$')

        with col3:
            options = ['exp1', 'exp2']
            datatable = st.selectbox('$\\bold{Select\ expression\ dataset:}$', options=options, index=None,
                                     placeholder='Choose a dataset',
                                     key='datasetwidget_%s' % species)

        if counts or datatable:
            col1, col2 = st.columns([0.4, 0.6])
            with col1:
                user_geneids = st.text_area('$\\bold{Introduce\ one\ or\ more\ geneIDs:}$',
                                            key='gidwidget_%s' % species)

            if user_geneids:
                st.markdown("<br></br>", unsafe_allow_html=True)

                gene_df = genelist_dataframe(counts, user_geneids.split())
                gene_df2 = gene_df.reset_index()  # pass replicates from index into a column
                gene_df_melt = pd.melt(gene_df2, id_vars=['Tissue', 'index'], value_vars=user_geneids.split())
                gene_df_melt.columns = ['Tissue', 'Replicates', 'Gene', 'Counts']

                with st.expander(label='$\Large Replicates$'):

                    col1, col2, col3 = st.columns([0.3, 0.4, 0.3])
                    with col2:
                        selected_gene = st.selectbox('Select gene:', options=user_geneids.split(), index=None,
                                                     key='glist_%s' % species)

                    if selected_gene:
                        st.markdown("<br></br>", unsafe_allow_html=True)
                        fig = px.scatter(gene_df_melt[gene_df_melt['Gene'] == selected_gene], x='Tissue', y='Counts',
                                         color='Replicates',
                                         hover_data='Replicates')
                        ##cbd6ce
                        fig.update_layout(plot_bgcolor="#495e5c",
                                          title={
                                              'text': "Expression values for {0} gene in different tissues".format(
                                                  selected_gene
                                              ),
                                              'y': 1.0,
                                              'x': 0.5,
                                              'xanchor': 'center',
                                              'yanchor': 'top'},
                                          yaxis_title="Counts",
                                          hovermode='x unified'
                                          )

                        fig.update_traces(marker_size=10, hoverlabel_bgcolor='blue')
                        st.plotly_chart(fig, use_container_width=True)
                    # "<h5 style='text-align: left; color: green;'> protein sequence(s):</h5>"

                with st.expander(label='$\Large Average\ values$'):
                    col1, col2, col3 = st.columns([0.3, 0.4, 0.3])
                    with col2:
                        selected_genes = st.multiselect('Select gene:', options=user_geneids.split(),
                                                        key='glistavg_%s' % species)

                    if selected_genes:
                        datafr = get_avg_value_per_tissue(gene_df_melt)

                        figr = px.line(datafr.loc[datafr['Gene'] == selected_genes[0]], x='Tissue', y='Average',
                                       color='Gene')
                        figr.update_traces(showlegend=False)
                        for i, n in enumerate(selected_genes):
                            temp_df = datafr.loc[datafr['Gene'] == n]
                            figr.add_scatter(
                                x=temp_df['Tissue'].values.tolist(),
                                y=temp_df['Average'].values.tolist(),
                                name=n
                            )

                        st.plotly_chart(figr, use_container_width=True)

                        # for g in selected_genes:
                        #     temp_df = datafr[datafr['Gene'] == g]
                        #     fig.add_trace(
                        #         go.Line(temp_df, x=temp_df['Tissue'], y=temp_df['Average'], color=temp_df['Gene'])
                        #     )

                with st.expander(label='$\Large Table$'):
                    st.markdown(
                        "<div style='text-align: center'> <b>Expression values table\n\n </b></div>",
                        unsafe_allow_html=True)

                    table_df = pd.DataFrame(gene_df[user_geneids.split()]).T
                    table_df = table_df.sort_index(axis=1)
                    st.write(table_df)
    with tab7:
        pass
    with tab8:
        pass
