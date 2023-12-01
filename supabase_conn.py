3#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Created on Tue Dec  3 08:26:27 2013

@author: kalexiou
"""

from supabase import create_client, Client
import streamlit as st
from io import StringIO

# Initialize connection.
url: str = st.secrets['connections']['supabase']["SUPABASE_URL"]
key: str = st.secrets['connections']['supabase']["SUPABASE_KEY"]
supabase: Client = create_client(url, key)

def getlist_from_fasta(uploaded_fasta):
    fasta_list = StringIO(uploaded_fasta.getvalue().decode("utf-8")).read().split('>')[1:]
    fasta_list2 = []
    for i in fasta_list:
        a = i.split('\n')
        tid = a[0].split(' ')[0]
        seq = ''.join(a[1:])
        fasta_list2.append([tid, seq])

    return fasta_list2


def import_data(datalist, species, dtype=('gene', 'cds', 'protein')):

    list_to_import = []
    for l in datalist:
        if dtype == 'gene':
            gid, chrom, start, end, transcripts, proteins, function = l
            list_to_import.append({"geneid": gid, "chrom": chrom, "start": start, "end": end,
                                   "transcripts": transcripts, "proteins":proteins, "function": function, "species": species})
        elif dtype == 'cds':
            recid, recseq = l
            list_to_import.append({"seqid": recid, "cdsFasta": str(recseq), "species": species})
        elif dtype == 'protein':
            recid, recseq = l
            list_to_import.append({"seqid": recid, "proteinFasta": str(recseq), "species": species})

    if dtype == 'gene':
        for g in list_to_import:
            supabase.table("gene_info").insert(g).execute()

    if dtype == 'cds':
        for g in list_to_import:
            supabase.table("cds_fastas").insert(g).execute()

    if dtype == 'protein':
        for g in list_to_import:
            supabase.table("protein_fastas").insert(g).execute()

    st.success('You imported your data!')

    return


# View Details
def view_all_data(tablename):
    # get all data from speciesdb table
    info, count = supabase.table(tablename).select("*").execute()
    all_data = list(info[1])

    data_list = []
    for d in all_data:
        data_list.append(list(d.values()))

    return data_list


def view_selected_data(tablename, sel_species, colname, searchterm):

    info, count = supabase.table(tablename).select('*').eq("species", sel_species).ilike('%s'%colname, '%{0}%'.format(searchterm)).execute()

    all_data = list(info[1])

    mlist = []
    for d in all_data:
        mlist.append(list(d.values())[1:])  # exclude index

    return mlist


def view_species_data(tablename, sel_species):
    info, count = supabase.table(tablename).select("*").eq("species", sel_species).execute()

    all_data = list(info[1])

    mlist = []
    for d in all_data:
        mlist.append(list(d.values())[1:])  # exclude index

    return mlist


class FileDownloader(object):
    """docstring for FileDownloader
    >>> download = FileDownloader(data,filename).download()

    """

    def __init__(self, data, filename="myfile"):
        super(FileDownloader, self).__init__()
        self.data = data
        self.filename = filename

    def download(self):
        b64 = base64.b64encode(self.data.encode()).decode()
        new_filename = "{}_{}.csv".format(self.filename, timestr)

        href = f'<a href="data:file/csv;base64,{b64}" download="{new_filename}">Click Here!!</a>'
        st.markdown(href, unsafe_allow_html=True)
