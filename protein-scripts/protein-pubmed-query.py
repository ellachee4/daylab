# Author: Ella Chee
# Date: July 2024

# Description: Script to query PubMed, given a CSV list of proteins, 
# returns CSV containing number of articles associated with each protein.
# uniprot id, check pubmed/pubmed central, terms within paper
# look at StringDB.org for protein-protein interactions
# look at antibody availability of proteins within paper
# look at protein size
# sort by range of years?

# Imports
from Bio import Entrez
import pandas as pd
import sys
import csv

Entrez.email = 'chee.el@northeastern.edu'

key_words = ['UV', 'UV damage', 'DNA damage', 'nucleotide excision repair', 'NER', 'DNA repair', 
             'G4', 'G-quadruplex', 'secondary structure', 'DNA-binding proteins', 'helicase']

def get_article_count(protein_name):
    '''
    Query PubMed for the number of articles associated with a given protein using both MeSH terms and text words.
    '''
    term = f'{protein_name}[MeSH Terms] OR {protein_name}[tw]'
    handle = Entrez.esearch(db='pubmed', term=term, keywds=key_words)
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
    
    # Sort list from smallest to largest
    protein_article_counts.sort(key=lambda x: x[1])

    return protein_article_counts

def main(input_csv, output_csv):
    '''
    Read a list of proteins from a CSV file, query PubMed for the 
    number of articles associated with each protein, save to CSV file.
    '''

    # Read input CSV
    prot_data = []
    with open(input_csv, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        for row in reader:
            prot_data.append(row)
    print('Read input CSV proteins')

    # Sorted list of proteins and their associated article counts
    sorted_protein_article_counts = query_pubmed(prot_data)
    print('Queried PubMed')

    # Save to CSV
    sorted_fields = ['Protein', 'Article Count']
    with open(output_csv, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(sorted_fields)
        for row in sorted_protein_article_counts:
            writer.writerow(row)
    print(f"Results saved to {output_csv}")

# List of proteins to query
proteins = ['POT1', 'RPA', 'CST', 'CTC1-STN1-TEN1', 'BRCA1', 'hnRNP A1', 'Pif1']

# Main: Adjust input and output CSV file names accordingly
if __name__ == "__main__":
    input_csv = "unique_to_mTbG4P_dedup.csv"
    output_csv = "protein_article_counts.csv"
    main(input_csv, output_csv)