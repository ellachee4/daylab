# Author: Ella Chee, Annika Salpukas
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
import requests
from collections import defaultdict

Entrez.email = 'salpukas.a@northeastern.edu'

key_words = ['UV', 'UV damage', 'DNA damage', 'nucleotide excision repair', 'NER', 'DNA repair',
             'G4', 'G-quadruplex', 'secondary structure', 'DNA-binding proteins', 'helicase']


def get_article_count(protein_name):
    '''
    Query PubMed for the number of articles associated with a given protein using both MeSH terms and text words.
    '''
    if protein_name=='NA':
        return 0
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
        protein_article_counts.append(count)
        length = len(protein_article_counts)
        if length % 100 == 0:
            print(length, '/', len(proteins), 'done')

    # Sort list from smallest to largest
    #protein_article_counts.sort(key=lambda x: x[1])
    print('Finished article count')

    return protein_article_counts


def uniprot_to_string(uniprots):
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
    map = {line.split("\t")[0]:line.split("\t")[2] for line in results.text.strip().split("\n")}

    # hard code mapping failures
    map['A8MXQ7'] = '9606.ENSP00000489685'
    map['P0CJ79'] = '9606.ENSP00000491567'
    map['Q8TD47'] = '9606.ENSP00000486252'

    print('Mapping successful')
    return map


def query_string(strings):
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
        l = line.strip().split("\t")
        str_dict[l[0]] += 1
        str_dict[l[1]] += 1

    print('Queried protein interactions')
    return str_dict


def go_score(go):
    term = go[1]
    positives = {'damage response', 'DNA binding', 'DNA-binding', 'DNA repair',
                 'nucleotide-excision repair', 'transcription'}
    negatives = {'RNA splicing', 'RNA processing', 'myosin', 'translation'}
    if any(word in term for word in positives):
        return 1
    elif any(word in term for word in negatives):
        return -1
    else:
        return 0


def query_go_terms(uniprots):
    scores = []
    terms = []
    for u in uniprots:
        url = f"https://www.ebi.ac.uk/proteins/api/proteins/{u}"
        headers = {
            "Accept": "application/json"
        }
        response = requests.get(url, headers=headers)
        data = response.json()
        gos = [(x['id'], x.get('properties', {}).get('term', 'N/A')) for x in data['dbReferences'] if x['type'] == 'GO']
        score = sum(map(go_score, gos))
        scores.append(score)
        terms.append([go[0] for go in gos])
    print('Queried GO terms')
    return scores, terms


def main(input_csv, output_csv):
    '''
    Read a list of proteins from a CSV file, query PubMed for the
    number of articles associated with each protein, save to CSV file.
    '''

    # Read input CSV
    uniprots = []
    symbols = []
    with open(input_csv, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        next(reader, None)
        for row in reader:
            uniprots.append(row[0])
            symbols.append(row[2])
    print('Read input CSV proteins')
    # Sorted list of proteins and their associated article counts
    article_counts = query_pubmed(symbols)
    map = uniprot_to_string(uniprots)
    print('Not converted:', [u for u in uniprots if u not in map.keys()])
    strings = [map[uni] if uni in map.keys() else '' for uni in uniprots]
    interactions = query_string(strings)
    go_scores, go_terms = query_go_terms(uniprots)

    # Save to CSV
    sorted_fields = ['Uniprot', 'Protein Symbol', 'Article Count', 'Interactions', 'GO Score', 'GO Terms']
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(sorted_fields)
        for i in range(len(article_counts)):
            writer.writerow([uniprots[i], symbols[i], article_counts[i],
                             interactions[strings[i]], go_scores[i], go_terms[i]])
    print(f"Results saved to {output_csv}")


# Main: Adjust input and output CSV file names accordingly
if __name__ == "__main__":
    input_csv = "unique_to_mTbG4P_dedup.csv"
    output_csv = "protein_article_counts2.csv"
    main(input_csv, output_csv)