# Author: Ella Chee, Annika Salpukas
# Date: July 2024

# Description: Script to query PubMed, given a CSV list of proteins,
# returns CSV containing number of articles associated with each protein 
# based on UNIPROT ID queried from PubMed, GO terms, protein-protein interactions from StringDB, 
# antibody availability, and final score for protein selection based on factors (weighted by importance).

# Imports
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
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
    else:
        term = f'{protein_name}[MeSH Terms] OR {protein_name}[tw] AND (UV[Title/Abstract] OR Ultraviolet radiation[Title/Abstract] OR G4[Title/Abstract] OR quadruplex[Title/Abstract] OR dna repair[Title/Abstract] OR melanoma[Title/Abstract])'
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
    '''Access StringDB API to convert Uniprot IDs to StringDB IDs.'''

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
    '''Access StringDB API to query protein-protein interactions.'''

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
    '''
    Calculate GO score, based on importance of GO term
    (extra weights for damage response, DNA binding, DNA repair, less for general repair).
    '''

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


#----------------- UV Score -----------------#
def map_uv(timepoint):
    '''Assign timepoints to UV scores. '''
    if timepoint=='mock':
        return -4
    elif timepoint=='1h':
        return 1
    elif timepoint=='24h':
        return 2
    else:
        return 3

def uv_score(uv):
    ''' Calculate UV score based on the timepoints of UV treatment.'''
    uv_list = list(set(uv.split(':')))
    print('UV scores calculated')
    return sum(map(map_uv, uv_list))

def query_go_terms(uniprots):
    '''Query GO terms for a list of proteins and return a list of GO scores and terms.'''

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
    '''
    Scrapes antibodypedia.com for the access link, number of antibodies, 
    number of referenced antibodies, and providers for a given UniProt ID.
    '''

    # Set up the Selenium  WebDriver, construct URL, navigate to URL
    options = Options()
    options.add_argument('--start-maximized')
    driver = webdriver.Chrome(options=options)
    url = f'https://www.antibodypedia.com/explore/uniprot%3A{uniprot_id}'
    try:
        driver.get(url)
    except Exception as error:
        print(f"Error for {uniprot_id}: {error}")

    # Load table
    WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.ID, "search_results_table"))
    )

    # Extract the link to the antibodies and number of antibodies
    try:
        link_tag = driver.find_element(By.XPATH, '//*[@id="search_results_table"]/tbody/tr/td[6]/a')
        if link_tag:
            antibodies_link = link_tag.get_attribute('href')
            antibody_id = str(link_tag.get_attribute('href'))[35:]
            number_of_antibodies = link_tag.text.split(' ')[0]
    except Exception as error:
        print('Error: no antibodies found for UNIPROT:', uniprot_id)
        antibodies_link = 'None'
        number_of_antibodies = 0
        antibody_id = 'None'
        
    # Extract number of providers from the text within the div
    try:
        providers_div = driver.find_element(By.CLASS_NAME, 'txtOne')
        if providers_div:
            try:
                number_of_providers = providers_div.find_element(By.XPATH, '//*[@id="search_results_table"]/tbody/tr/td[6]/div/b').text
            except Exception as error:
                number_of_providers = 0
    except Exception as error:
        number_of_providers = 0

    driver.quit()
    print('Queried for protein:', uniprot_id)
    return (antibodies_link, antibody_id, number_of_antibodies, number_of_providers)

def track_references_antipodypedia(antibody_id):
    '''
    Track antibodies with references from antibodypedia.
    '''

    if antibody_id == 'None':
        return 0
    
    # Set up the Selenium  WebDriver, construct URL, navigate to URL
    driver = webdriver.Chrome()
    url = f'https://www.antibodypedia.com/gene/{antibody_id}?reference%5B%5D=yes'
    driver.get(url)
    
    # Load table
    WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.ID, "content"))
    )

    # Extract the link to the antibodies and number of antibodies
    try: 
        link_tag = driver.find_element(By.XPATH, '//*[@id="filter_results"]/b[1]')
        if link_tag:
            referenced_antibodies = link_tag.text
    except Exception as error:
        print('Error: no references found for antibody:', antibody_id)
        referenced_antibodies = 0
    driver.quit()
    return referenced_antibodies

def query_antibodypedia(uniprots):
    '''
    Query antibodypedia given a list of proteins and return a list of 
    antibody links, number of antibodies (including referenced), and number of providers.
    '''

    links = list()
    ids = list()
    antibodies = list()
    providers = list()
    references = list()

    for u in uniprots:
        link, id, antibody, provider = scrape_antibodypedia_data(u)
        links.append(link)
        ids.append(id)
        antibodies.append(antibody)
        providers.append(provider)
    
    for i in ids:
        references.append(track_references_antipodypedia(i))
    print('Successfully queried antibodypedia')
    return links, antibodies, providers, references

#----------------- Compute Score -----------------#
def compute_score(row):
    '''
    Compute weighted score for protein selection based on article count (favor fewer), 
    number of referenced antibodies, interactions GO score, and UV score.
    '''

    # late more important, mock less important ercc6 rated high by score, 
    # present at 3 timepoints, mock is one of them (fix uv score)
    # any term present gets increase score if found in abstract

    score = row['Article Count (normalized)']*-2 + row['Number of ref-ed Antibodies (normalized)']*0.5\
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
    links, num_antibodies, num_providers, num_referenced = query_antibodypedia(uniprots)
    article_counts = query_pubmed(symbols)
    mapping = uniprot_to_string(uniprots)
    print('Not converted:', [u for u in uniprots if u not in mapping.keys()])
    strings = [mapping[uni] if uni in mapping.keys() else '' for uni in uniprots]
    interactions = query_string(strings)
    go_scores, go_terms = query_go_terms(uniprots)

    data_dict = {
        'Uniprot': uniprots,
        'Protein Symbol': symbols,
        'Article Count': article_counts,
        'Interactions': [interactions[string] for string in strings],
        'Antibody Link': links,
        'Number of Antibodies': num_antibodies,
        'Number of Providers': num_providers,
        'Number of Ref-ed Antibodies': num_referenced,
        'UV_treatment': uvs,
        'GO Score': go_scores,
        'GO Terms': go_terms
    }

    df_all = pd.DataFrame(data_dict)
    df_all['UV Score'] = df_all['UV_treatment'].apply(lambda x: uv_score(x))
    for col in ['Number of Ref-ed Antibodies', 'Article Count', 'Interactions', 'GO Score', 'UV Score']:
        df_all[col] =pd.to_numeric(df_all[col])
        df_all[f'{col} (normalized)'] = df_all[col]/df_all[col].std()

    # Compute weighted score for protein selection and save results to CSV
    df_all['Overall Score'] = df_all.apply(lambda x: compute_score(x), axis=1)
    cols = ['Uniprot', 'Protein Symbol','Antibody Link', 'Number of Antibodies', 'Number of Providers', 'Number of Ref-ed Antibodies',
            'Article Count', 'Interactions', 'GO Score', 'UV Score','Overall Score', 'GO Terms']
    df_all = df_all[cols]
    df_all.to_csv(final_output_csv)
    print('Protein selection data saved to', final_output_csv)

# Main: Adjust input and output CSV file names accordingly
if __name__ == "__main__":
    input_csv = "proteins_with_uv.csv"
    final_output_csv = "combined-scores3.csv"
    main(input_csv, final_output_csv)