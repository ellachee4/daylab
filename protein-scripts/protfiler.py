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


def get_article_count(protein_name, do_all):
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
        term = f'{protein_name}[MeSH Terms] OR {protein_name}[tw] AND (UV[Title/Abstract] OR Ultraviolet radiation[Title/Abstract] OR G4[Title/Abstract] OR quadruplex[Title/Abstract] OR dna repair[Title/Abstract] OR melanoma[Title/Abstract])'

    handle = Entrez.esearch(db='pubmed', term=term)
    record = Entrez.read(handle)
    handle.close()
    return int(record['Count'])


def query_pubmed(proteins):
    """
    Query PubMed for a list of proteins and return a sorted list of proteins
    and the number of associated articles.
    """
    protein_all_article_counts = []
    protein_specific_article_counts = []
    for protein in proteins:
        all_count = get_article_count(protein, True)
        protein_all_article_counts.append(all_count)
        specific_count = get_article_count(protein, False)
        protein_specific_article_counts.append(specific_count)
        length = len(protein_all_article_counts)
        if length % 100 == 0:
            print(length, '/', len(proteins), 'done')

    print('Finished article count')
    return protein_all_article_counts, protein_specific_article_counts


# ----------------- Protein-Protein Interactions ----------------- #


def uniprot_to_string(uniprots):
    """Access StringDB API to convert Uniprot IDs to StringDB IDs."""

    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "get_string_ids"

    params = {
        "identifiers": "\r".join(uniprots),
        "species": 9606,
        "limit": 1,
        "echo_query": 1
    }

    request_url = "/".join([string_api_url, output_format, method])
    results = requests.post(request_url, data=params)
    mapper = {line.split("\t")[0]: line.split("\t")[2] for line in results.text.strip().split("\n")}

    # Hard code mapping failures
    mapper['A8MXQ7'] = '9606.ENSP00000489685'
    mapper['P0CJ79'] = '9606.ENSP00000491567'
    mapper['Q8TD47'] = '9606.ENSP00000486252'

    print('Mapping successful')
    return mapper


def query_string(strings):
    """Access StringDB API to query protein-protein interactions."""

    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"

    params = {
        "identifiers": "%0d".join(strings),
        "species": 9606
    }

    request_url = "/".join([string_api_url, output_format, method])
    response = requests.post(request_url, data=params)

    lines = set(response.text.strip().split("\n"))
    str_dict = defaultdict(lambda: 0)
    for line in lines:
        line_list = line.strip().split("\t")
        if len(line_list) >= 2:
            str_dict[line_list[0]] += 1
            str_dict[line_list[1]] += 1

    print('Queried protein interactions')
    return str_dict


# ----------------- GO Terms ----------------- #


def query_ebi(uniprots):
    """Query GO terms for a list of proteins and return a list of GO scores and terms."""

    scores = []
    terms = []
    antibody_ids = []
    for u in uniprots:

        # Define url for Protein API
        url = f"https://www.ebi.ac.uk/proteins/api/proteins/{u}"
        headers = {
            "Accept": "application/json"
        }
        response = requests.get(url, headers=headers)
        data = response.json()
        gos = [(x['id'], x.get('properties', {}).get('term', 'N/A')) for x in data['dbReferences'] if x['type'] == 'GO']
        antibodypedia_id = ''

        # Collect GO terms and Antibodypedia IDs
        for x in data['dbReferences']:
            if x['type'] == 'GO':
                gos.append((x['id'], x.get('properties', {}).get('term', 'N/A')))
            elif x['type'] == 'Antibodypedia':
                antibodypedia_id = x['id']
        score = sum(map(go_score, gos))
        scores.append(score)
        terms.append([go[1] for go in gos])
        antibody_ids.append(antibodypedia_id)

    print('Queried GO terms')
    return scores, terms, antibody_ids


# ----------------- Antibodypedia ----------------- #


def scrape_antibodypedia(symbol, antibody_id):
    """
    Scrapes Antibodypedia to return number of referenced antibodies
    for a particular protein, return link and number
    """

    if symbol != 'NA' and antibody_id != '':

        # Define URL and create BeautifulSoup object
        url = f"https://www.antibodypedia.com/gene/{antibody_id}/{symbol}?reference%5B%5D=yes"
        response = requests.get(url)
        soup = BeautifulSoup(response.text, 'html.parser')

        div = soup.find(id='filter_results')
        if div:
            bold_tags = div.find_all('b')

            if bold_tags:
                return f"https://www.antibodypedia.com/gene/{antibody_id}/{symbol}", int(
                    bold_tags[0].get_text(strip=True))
            else:
                return 'None', 0
        else:
            return 'None', 0
    else:
        return 'None', 0


def query_antibodypedia(symbols, antibody_ids):
    """
    Query Antibodypedia given a list of proteins and return a list of
    antibody links, and number of referenced antibodies.
    """

    links = list()
    num_referenced = list()

    # Iterate through symbols and ids
    for sym, abid in zip(symbols, antibody_ids):
        link, n_ref = scrape_antibodypedia(sym, abid)
        links.append(link)
        if n_ref > 2:
            ref_score = 2
        elif n_ref > 0:
            ref_score = 1
        else:
            ref_score = 0
        num_referenced.append(ref_score)
    print('Queried Antibodypedia')
    return links, num_referenced


# ----------------- GO Score ----------------- #


def go_score(go):
    """
    Calculate GO score, based on importance of GO term
    (extra weights for damage response, DNA binding, DNA repair, less for general repair).
    """

    # Define positive and negative terms
    term = go[1]
    positives = {'damage response', 'DNA binding', 'DNA-binding', 'DNA repair',
                 'nucleotide-excision repair', 'transcription', 'repair', 'DNA helicase', 'helicase',
                 'chromatin binding', 'G4', 'quadruplex', 'guanine', 'gtpase', 'transcription initiation',
                 'transcription termination', 'transcription activator', 'poly-ADP-D-ribose', 'nuclear',
                 'ubiquitin', 'melanoma', 'autophagy', 'apoptosis', 'replication', 'damaged DNA binding',
                 'nucleosome', 'histone', 'regulatory'}

    negatives = {'RNA splicing', 'RNA processing', 'myosin', 'translation', 'ribosomal', 'ribosome',
                 'cytosol', 'cytosolic', 'keratinization -2'}

    # Return scores
    if any(word in term for word in positives):
        return 1
    elif any(word in term for word in negatives):
        return -1
    else:
        return 0

# ----------------- UV Score ----------------- #


def map_uv(timepoint):
    """Assign timepoints to UV scores. """
    if timepoint == 'mock':
        return -4
    elif timepoint == '01h':
        return 1
    elif timepoint == '24h':
        return 2
    else:
        return 3


def uv_score(uv):
    """ Calculate UV score based on the timepoints of UV treatment."""
    uv_list = list(set(uv.split(':')))
    return sum(map(map_uv, uv_list))


# ----------------- Compute Score ----------------- #


def compute_score(row, cols, weights):
    """
    Compute weighted score for protein selection based on article counts (favor fewer),
    number of referenced antibodies, interactions GO score, and UV score.
    """
    score = sum(row[col + ' (normalized)'] * weight for col, weight in zip(cols, weights))
    return score


# ----------------- Main ----------------- #


def run_protfiler(args):
    """
    Read a list of proteins from a CSV file, query PubMed for the
    number of articles associated with each protein, string interactions,
    query Antibodypedia.com for referenced antibodies, find associated
    GO terms and compute a selection score, save results to CSV file.
    """

    # Track start time
    start_time = time.time()

    # Initialize lists
    uniprots = []
    symbols = []
    uvs = []

    # Get indices for id columns
    uniprot_index = int(args.uni)
    symbol_index = int(args.sym)
    uv_index = int(args.uv) if args.uv else None

    # Open and read input CSV
    with open(args.i, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        next(reader, None)
        for row in reader:
            uniprots.append(row[uniprot_index])
            symbols.append(row[symbol_index])
            if args.uv:
                uvs.append(row[uv_index])
    print('Read input CSV proteins')

    # Sorted list of proteins and their associated article counts
    total_article_counts, term_article_counts = query_pubmed(symbols)
    mapping = uniprot_to_string(uniprots)
    print('String IDs not found:', [u for u in uniprots if u not in mapping.keys()])
    strings = [mapping[uni] if uni in mapping.keys() else '' for uni in uniprots]
    interactions = query_string(strings)
    go_scores, go_terms, antibodypedia_ids = query_ebi(uniprots)
    links, num_referenced = query_antibodypedia(symbols, antibodypedia_ids)

    # Initialize data dictionary
    data_dict = {
        'Uniprot': uniprots,
        'Protein Symbol': symbols,
        'Total Article Count': total_article_counts,
        'Term Article Count': term_article_counts,
        '\'String\' Interactions': [interactions[string] for string in strings],
        'Antibodypedia Link': links,
        'Ref-ed Antibodies via Antibodypedia': num_referenced,
        'GO Score': go_scores,
        'GO Terms': go_terms
    }

    # Set up dataframe and normalize counts/scores
    df_all = pd.DataFrame(data_dict)
    score_cols = ['Ref-ed Antibodies via Antibodypedia', 'Total Article Count', 'Term Article Count',
                  '\'String\' Interactions', 'GO Score']
    score_weights = [1, -1, 1, 1.5, 1]

    # In case of UV timepoints
    if args.uv:
        df_all['UV_treatment'] = uvs
        df_all['UV Score'] = df_all['UV_treatment'].apply(lambda x: uv_score(x))
        score_cols.append('UV Score')
        score_weights.append(2)

    # Normalize columns for scoring
    for col in score_cols:
        df_all[col] = pd.to_numeric(df_all[col])
        df_all[f'{col} (normalized)'] = df_all[col]/df_all[col].std()

    # Compute weighted score for protein selection
    df_all['Overall Score'] = df_all.apply(lambda x: compute_score(x, score_cols, score_weights), axis=1)
    print('Computed scores')
    cols = ['Uniprot', 'Protein Symbol', 'Antibodypedia Link'] + score_cols + ['Overall Score', 'GO Terms']
    df_all = df_all[cols]

    # Sort by overall score and save to CSV
    df_all.sort_values(by=['Overall Score'], ascending=False, inplace=True)
    df_all.to_csv(args.o)

    # Track end time
    endtime = time.time()
    total_time = endtime - start_time
    print(f'Time taken: {total_time:.2f} seconds')
    print('Protein selection data saved to', args.o)


def main():
    """
    Main function reads command line arguments and runs protfiler
    """

    # Parse command line arguments
    parser = ArgumentParser(description='Run Protfiler')
    parser.add_argument('-i', help='Input csv file')
    parser.add_argument('-o', help='Output csv file')
    parser.add_argument('-uni', help='Input uniprot column index')
    parser.add_argument('-sym', help='Input symbol column index')
    parser.add_argument('-uv', help='Input uv timepoints column index')

    args = parser.parse_args()

    # Run Protfiler
    run_protfiler(args)


# Main: Adjust input and output CSV file names accordingly on command line
if __name__ == "__main__":
    main()
