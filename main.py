# This is a sample Python script.

# Press Mayús+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import pandas as pd
import base64
from page_layout import *

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

        if selected_species == 'cannabis :herb:':
            if 'keyword' not in st.session_state or 'idsearch' not in st.session_state:
                st.session_state.keyword = ''
                st.session_state.idsearch = ''
            generate_page(selected_species)

        elif selected_species == 'b':
            if 'keyword' not in st.session_state or 'idsearch' not in st.session_state:
                st.session_state.keyword = ''
                st.session_state.idsearch = ''
            generate_page(selected_species)


if __name__ == "__main__":
    main()
