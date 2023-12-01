####2023-11-16
1. Obtaining cannabis gene information.
 
    - select genes

                zcat GCF_900626175.2_cs10_genomic.gff.gz | awk '$3=="gene"' | cut -f9 | cut -d';' -f1 > geneids  # 29,813

####2023-11-17
1. Extract gff entries for genes and the associated RNA information.

            zgrep -v '#' GCF_900626175.2_cs10_genomic.gff.gz | awk '$3!="exon" && $3!="CDS"' > intermed.gff

            while read line; do echo $line; grep -A1 $line intermed.gff; done <geneids >> geneids_plus-rnas.tab

####2023-11-20
1. New strategy applied: I have written a python script for obtaininging gene coordinates, together with rnas and functions.

            python3 extractinfo.py

2. I moved into pycharm and started preparing the interface's back- and front-end.
