# Author: Ella Chee, Annika Salpukas
# Date: July 2024

# Description: Script to query PubMed, given a CSV list of proteins,
# returns CSV containing number of articles associated with each protein 
# based on UNIPROT ID queried from PubMed, GO terms, protein-protein interactions from StringDB, 
# antibody availability, and final score for protein selection based on factors (weighted by importance).

# Imports
import sys
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from Bio import Entrez
import sys
import csv
import requests
from collections import defaultdict
import ast

#----------------- Article Count -----------------#

# Adjust email as needed
Entrez.email = 'salpukas.a@northeastern.edu'

def get_article_count(protein_name):
    '''
    Query PubMed for the number of articles associated with a 
    given protein using both MeSH terms and text words.
    '''

    if protein_name=='NA':
        return 0
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
        protein_article_counts.append(count)
        length = len(protein_article_counts)
        if length % 100 == 0:
            print(length, '/', len(proteins), 'done')

    print('Finished article count')
    return protein_article_counts

#----------------- Protein-Protein Interactions -----------------#
def uniprot_to_string(uniprots):
    '''Access StringDB API to convert Uniprot IDs to StringDB IDs'''

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
    '''Access StringDB API to query protein-protein interactions'''

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

#----------------- GO Terms -----------------#

def go_score(go):
    '''Calculate GO score, based on importance of GO term
    (extra weights for damage response, DNA binding, DNA repair, less for general repair)'''

    term = go[1]
    positives = {'damage response', 'DNA binding', 'DNA-binding', 'DNA repair',
                 'nucleotide-excision repair', 'transcription', 'repair'}
    negatives = {'RNA splicing', 'RNA processing', 'myosin', 'translation'}
    # extra weight for damage response, dna binding, dna repair, less weight for general repair
    if any(word in term for word in positives):
        return 1
    elif any(word in term for word in negatives):
        return -1
    else:
        return 0

def uv_score(uv):
    uv_list = list(set(uv.split(':')))
    return sum(map(lambda x: -1 if x=='mock' else 1, uv_list))

def query_go_terms(uniprots):
    '''Query GO terms for a list of proteins and return a list of GO scores and terms'''

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
        terms.append([go[1] for go in gos])
    print('Queried GO terms')
    return scores, terms

#----------------- Antibody Availability -----------------#
def scrape_antibodypedia_data(uniprot_id):
    '''Scrapes antibodypedia.com for the access link, 
    number of antibodies and providers for a given UniProt ID'''

    # Set up the Selenium  WebDriver, construct URL, navigate to URL
    driver = webdriver.Chrome()
    base_url = 'https://www.antibodypedia.com/explore/uniprot%3A'
    url = f'{base_url}{uniprot_id}'
    driver.get(url)

    # Load table
    WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.ID, "search_results_table"))
    )
    print('Found table for UNIPROT:', uniprot_id)

    # Extract the link to the antibodies and number of antibodies
    try:
        link_tag = driver.find_element(By.XPATH, '//*[@id="search_results_table"]/tbody/tr/td[6]/a')
        if link_tag:
            antibodies_link = link_tag.get_attribute('href')
            number_of_antibodies = link_tag.text.split(' ')[0]
            print('Antibodies Link:', antibodies_link)
            print('Number of Antibodies:', number_of_antibodies)
    except Exception as error:
        print('Error: no antibodies found for UNIPROT:', uniprot_id)
        antibodies_link = 'None'
        number_of_antibodies = 0
        
    # Extract number of providers from the text within the div
    try:
        providers_div = driver.find_element(By.CLASS_NAME, 'txtOne')
        if providers_div:
            try:
                number_of_providers = providers_div.find_element(By.XPATH, '//*[@id="search_results_table"]/tbody/tr/td[6]/div/b').text
                print('Number of Providers:', number_of_providers)
            except Exception as error:
                print('Error: no providers found for UNIPROT:', uniprot_id)
                number_of_providers = 0
    except Exception as error:
        print('Error: no providers found for UNIPROT:', uniprot_id)
        number_of_providers = 0


    driver.quit()
    return (antibodies_link, number_of_antibodies, number_of_providers)

def query_antibodypedia(uniprots):
    '''Query antibodypedia given a list of proteins and return 
    a list of antibody links, number of antibodies, and number of providers'''
    links = list()
    antibodies = list()
    providers = list()

    for u in uniprots:
        link, antibody, provider = scrape_antibodypedia_data(u)
        links.append(link)
        antibodies.append(antibody)
        providers.append(provider)

    return links, antibodies, providers

def compute_score(row):
    score = row['Article Count (normalized)']*-2 + row['Number of Antibodies (normalized)']*0.5\
            + row['Interactions (normalized)'] + row['GO Score (normalized)'] + row['UV Score (normalized)']
    return score

#----------------- Main -----------------#
def main(input_csv, output_csv):
    '''
    Read a list of proteins from a CSV file, query PubMed for the
    number of articles associated with each protein, save to CSV file.
    '''

    # Read input CSV
    uniprots = []
    symbols = []
    uvs = []
    with open(input_csv, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        next(reader, None)
        for row in reader:
            uniprots.append(row[0])
            symbols.append(row[2])
            uvs.append(row[10])
    print('Read input CSV proteins')

    # Sorted list of proteins and their associated article counts
    article_counts = query_pubmed(symbols)
    mapping = uniprot_to_string(uniprots)
    print('Not converted:', [u for u in uniprots if u not in mapping.keys()])
    strings = [mapping[uni] if uni in mapping.keys() else '' for uni in uniprots]
    interactions = query_string(strings)
    go_scores, go_terms = query_go_terms(uniprots)
    links, num_antibodies, num_providers = query_antibodypedia(uniprots)

    data_dict = {
        'Uniprot': uniprots,
        'Protein Symbol': symbols,
        'Article Count': article_counts,
        'Interactions': [interactions[string] for string in strings],
        'Antibody Link': links,
        'Number of Antibodies': num_antibodies,
        'Number of Providers': num_providers,
        'UV_treatment': uvs,
        'GO Score': go_scores,
        'GO Terms': go_terms
    }

    df_all = pd.DataFrame(data_dict)

    '''
    sorted_fields = ['Uniprot', 'Protein Symbol', 'Article Count', 'Interactions', 'GO Score', 'GO Terms']
    with open(go_output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(sorted_fields)
        for i in range(len(article_counts)):
            writer.writerow([uniprots[i], symbols[i], article_counts[i],
                             interactions[strings[i]], go_scores[i], go_terms[i]])
    '''
    df_all['UV Score'] = df_all['UV_treatment'].apply(lambda x: uv_score(x))
    for col in ['Number of Antibodies', 'Article Count', 'Interactions', 'GO Score', 'UV Score']:
        df_all[col] =pd.to_numeric(df_all[col])
        df_all[f'{col} (normalized)'] = df_all[col]/df_all[col].std()

    # Compute weighted score for protein selection

    df_all['Overall Score'] = df_all.apply(lambda x: compute_score(x), axis=1)
    cols = ['Uniprot', 'Protein Symbol','Antibody Link', 'Number of Antibodies', 'Number of Providers',
            'Article Count', 'Interactions', 'GO Score', 'UV Score','Overall Score', 'GO Terms']

    df_all = df_all[cols]

    df_all.to_csv(final_output_csv)
    print('Protein selection data saved to', final_output_csv)

# Main: Adjust input and output CSV file names accordingly
if __name__ == "__main__":
    input_csv = "proteins_with_uv.csv"
    final_output_csv = "combined-scores2.csv"
    main(input_csv, final_output_csv)