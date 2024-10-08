# Author: Ella Chee, Annika Salpukas
# Date: August 2024

# Description: Script to query PubMed, given a CSV list of proteins,
# returns CSV containing number of articles associated with each protein
# based on UNIPROT ID queried from PubMed, GO terms, protein-protein interactions from StringDB,
# antibody availability, and final score for protein selection based on factors (weighted by importance).

# import libraries
from argparse import ArgumentParser
import pandas as pd
import time
from Bio import Entrez
import csv
from collections import defaultdict
import requests
from bs4 import BeautifulSoup

# Adjust email as needed
Entrez.email = 'salpukas.a@northeastern.edu'

# ----------------- Article Count ----------------- #

def get_article_count(protein_name, article_terms, do_all):
    """
    Query PubMed for the number of articles associated with a
    given protein using both MeSH terms and text words.
    """
    if protein_name == 'NA':
        return 0

    # Query for articles associated with the protein name (do_all=True) or
    # articles associated with the protein name and specific terms (do_all=False)
    if do_all:
        term = f'{protein_name}[MeSH Terms] OR {protein_name}[tw]'
    else:
        joined = ' OR '.join([t + '[Title/Abstract]' for t in article_terms])
        term = f'{protein_name}[MeSH Terms] OR {protein_name}[tw] AND ({joined})'

    handle = Entrez.esearch(db='pubmed', term=term)
    record = Entrez.read(handle)
    handle.close()
    return int(record['Count'])


def query_pubmed(proteins, terms):
    """
    Query PubMed for a list of proteins and return a sorted list of proteins
    and the number of associated articles.
    """
    protein_all_article_counts = []
    protein_specific_article_counts = []
    for protein in proteins:
        all_count = get_article_count(protein, terms, True)
        protein_all_article_counts.append(all_count)
        specific_count = get_article_count(protein, terms, False)
        protein_specific_article_counts.append(specific_count)
        length = len(protein_all_article_counts)
        if length % 50 == 0:
            time.sleep(3)
        if length % 100 == 0:
            print(length, '/', len(proteins), 'done')

    print('Finished article count')
    return protein_all_article_counts, protein_specific_article_counts


# ----------------- GO Terms ----------------- #



# ----------------- GO Score ----------------- #


def go_score(go, positives, negatives):
    """
    Calculate GO score, based on importance of GO term
    (extra weights for damage response, DNA binding, DNA repair, less for general repair).
    """

    # Define positive and negative terms
    term = go[1]

    # Return scores
    if any(word in term for word in positives):
        return 1
    elif any(word in term for word in negatives):
        return -1
    else:
        return 0


# ----------------- Compute Score ----------------- #


def compute_score(row, cols, weights):
    """
    Compute weighted score for protein selection based on article counts (favor fewer),
    number of referenced antibodies, interactions GO score, and UV score.
    """
    score = sum(row[col + ' (normalized)'] * weight for col, weight in zip(cols, weights))
    return score


# ----------------- Main ----------------- #

def get_go_terms():
    """
    Retrieve positively and negatively weighted GO keywords from user
    """
    positives = input('Enter positively weighted GO keywords (comma-separated): ')
    positives = [word.strip() for word in positives.split(',')]
    negatives = input('Enter negatively weighted GO keywords (comma-separated): ')
    negatives = [word.strip() for word in negatives.split(',')]
    return positives, negatives


# CHANGE 
def get_terms_and_weights(uv):
    """
    Retrieve relevant article search terms and score weights from user
    """
    terms = input('Enter relevant terms for article search (comma-separated): ')
    terms = [word.strip() for word in terms.split(',')]

    default_weights = input('Use default weights? (Yes/No): ')
    weights = [1, -1, 1, 1.5, 1]
    question = """
    Enter weights for antibody availability, total article count, 
    relevant term article count, protein-protein interactions, 
    """
    if uv:
        weights.append(2)
        question += 'GO term score, and UV score (comma-separated): '
    else:
        question += 'and GO term score (comma-separated)'

    if default_weights != 'Yes':
        weights = input(question)
        weights = [float(w.strip()) for w in weights.split(',')]
    return terms, weights


def run_protfiler(args):
    """
    Read a list of proteins from a CSV file, query PubMed for the
    number of articles associated with each protein, string interactions,
    query Antibodypedia.com for referenced antibodies, find associated
    GO terms and compute a selection score, save results to CSV file.
    """
    positive_gos, negative_gos = get_go_terms()
    article_terms, score_weights = get_terms_and_weights(args.uv)

    # Track start time
    start_time = time.time()

    # Initialize lists
    symbols = []

    # Get indices for id columns
    symbol_index = int(args.sym)

    # Open and read input CSV
    with open(args.i, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        next(reader, None)
        for row in reader:
            symbols.append(row[symbol_index])
    print('Read input CSV proteins')

    # Sorted list of proteins and their associated article counts
    total_article_counts, term_article_counts = query_pubmed(symbols, article_terms)

    # FIX GO TERMS
    go_scores, go_terms, antibodypedia_ids = query_ebi(uniprots, positive_gos, negative_gos)

    # Initialize data dictionary
    data_dict = {
        'Gene Symbol': symbols,
        'Total Article Count': total_article_counts,
        'Term Article Count': term_article_counts,
        'GO Score': go_scores,
        'GO Terms': go_terms
    }

    # Set up dataframe and normalize counts/scores
    df_all = pd.DataFrame(data_dict)
    score_cols = ['Total Article Count', 'Term Article Count','GO Score']


    # Normalize columns for scoring
    for col in score_cols:
        df_all[col] = pd.to_numeric(df_all[col])
        if df_all[col].nunique() == 1:
            df_all[f'{col} (normalized)'] = df_all[col]
        else:
            df_all[f'{col} (normalized)'] = df_all[col]/df_all[col].std()

    # Compute weighted score for protein selection
    df_all['Overall Score'] = df_all.apply(lambda x: compute_score(x, score_cols, score_weights), axis=1)
    print('Computed scores')
    cols = ['Gene Symbol'] + score_cols + ['Overall Score', 'GO Terms']
    df_all = df_all[cols]

    # Sort by overall score and save to CSV
    df_all.sort_values(by=['Overall Score'], ascending=False, inplace=True)
    df_all.to_csv(args.o)

    # Track end time
    endtime = time.time()
    total_time = endtime - start_time
    print(f'Time taken: {total_time:.2f} seconds')
    print('Gene selection data saved to', args.o)


def main():
    """
    Main function reads command line arguments and runs protfiler
    """

    # Parse command line arguments
    parser = ArgumentParser(description='Run Protfiler')
    parser.add_argument('-i', help='Input csv file')
    parser.add_argument('-o', help='Output csv file')
    parser.add_argument('-sym', help='Input symbol column index')

    args = parser.parse_args()

    # Run Protfiler
    run_protfiler(args)

# Main: Adjust input and output CSV file names accordingly on command line
if __name__ == "__main__":
    main()
