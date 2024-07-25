from Bio import Entrez
import pandas as pd

Entrez.email = 'chee.el@northeastern.edu'

def get_article_count(protein_name):
    '''
    Query PubMed for the number of articles associated with a given protein using both MeSH terms and text words.
    '''
    term = f'{protein_name}[MeSH Terms] OR {protein_name}[tw]'
    handle = Entrez.esearch(db='pubmed', term=term)
    record = Entrez.read(handle)
    handle.close()
    return int(record['Count'])

def query_pubmed(proteins):
    '''
    Query PubMed for a list of proteins and return a sorted list of proteins 
    and the number of associated articles.
    '''
    protein_article_counts = []
    for protein in proteins:
        count = get_article_count(protein)
        protein_article_counts.append((protein, count))
    
    # sort list smallest to largest
    protein_article_counts.sort(key=lambda x: x[1])
    
    return protein_article_counts

def main(list_proteins):
    # sorted list of proteins and their associated article counts
    sorted_protein_article_counts = query_pubmed(list_proteins)

    # convert to dataframe
    df = pd.DataFrame(sorted_protein_article_counts, columns=['Protein', 'Article Count'])
    print(df)

# list of proteins to query
proteins = ['POT1', 'RPA', 'CST', 'CTC1-STN1-TEN1', 'BRCA1', 'hnRNP A1', 'Pif1']

main(proteins)