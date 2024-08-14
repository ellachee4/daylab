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
from selenium.webdriver.chrome.options import Options
from selenium.common.exceptions import TimeoutException
import time
from Bio import Entrez
import sys
import csv
import requests
from collections import defaultdict
import ast

#----------------- Article Count -----------------#

# Adjust email as needed
Entrez.email = 'chee.el@northeastern.edu'

def get_article_count(protein_name, doAll):
    '''
    Query PubMed for the number of articles associated with a 
    given protein using both MeSH terms and text words.
    '''
    if protein_name=='NA':
        return 0
    
    # Query for articles associated with the protein name (doAll=True) or
    # articles associated with the protein name and specific terms (doAll=False)
    if doAll:
        term = f'{protein_name}[MeSH Terms] OR {protein_name}[tw]'
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

    # Hard code mapping failures
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
        if len(l) >= 2:
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
                 'nucleotide-excision repair', 'transcription', 'repair', 'DNA helicase', 'helicase', 
                 'chromatin binding', 'G4', 'quadruplex', 'guanine', 'gtpase', 'transcription initiation',
                 'transcription termination', 'transcription activator', 'poly-ADP-D-ribose', 'nuclear',
                 'ubiquitin', 'melanoma', 'autophagy', 'apoptosis', 'replication', 'damaged DNA binding',
                 'nucleosome', 'histone', 'regulatory'}

    negatives = {'RNA splicing', 'RNA processing', 'myosin', 'translation', 'ribosomal', 'ribosome',
                 'cytosol', 'cytosolic', 'keratinization -2'}

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
    elif timepoint=='01h':
        return 1
    elif timepoint=='24h':
        return 2
    else:
        return 3

def uv_score(uv):
    ''' Calculate UV score based on the timepoints of UV treatment.'''
    uv_list = list(set(uv.split(':')))
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
    Scrapes antibodypedia.com for the access link, 
    and number of referenced antibodies, for a given UniProt ID.
    '''

    # Set up the Selenium  WebDriver, construct URL, navigate to URL
    options = Options()
    options.add_argument("--headless")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--dns-prefetch-disable")
    driver = webdriver.Chrome(options=options)
    url = f'https://www.antibodypedia.com/explore/uniprot%3A{uniprot_id}'
    try:
        driver.get(url)
    except Exception as error:
        driver.navigate().refresh()
        print(f"Error for {uniprot_id}: {error}")

    # Load table
    WebDriverWait(driver, 20).until(
        EC.presence_of_element_located((By.ID, "search_results_table"))
    )

    # Extract the link to the antibodies and number of antibodies
    try:
        link_tag = driver.find_element(By.XPATH, '//*[@id="search_results_table"]/tbody/tr/td[6]/a')
        if link_tag:
            antibodies_link = link_tag.get_attribute('href')
            antibody_id = str(link_tag.get_attribute('href'))[35:]
    except Exception as error:
        print('Error: no antibodies found for UNIPROT:', uniprot_id)
        antibody_id = 'None'
        antibodies_link = 'None'

    driver.quit()
    return (antibodies_link, antibody_id)

def track_references_antipodypedia(antibody_id):
    '''
    Track number antibodies with references from antibodypedia.
    '''

    # Check if there is no antibody found
    if antibody_id == 'None':
        return 0
    
    # Set up the Selenium  WebDriver, construct URL, navigate to URL
    options = Options()
    options.add_argument("--headless")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--dns-prefetch-disable")
    driver = webdriver.Chrome(options=options)
    url = f'https://www.antibodypedia.com/gene/{antibody_id}?reference%5B%5D=yes'
    driver.get(url)
    
    # Load table
    WebDriverWait(driver, 20).until(
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
    antibody links, and number of referenced antibodies.
    '''

    links = list()
    ids = list()
    num_referenced = list()

    for u in uniprots:
        link, id = scrape_antibodypedia_data(u)
        links.append(link)
        ids.append(id)
    
    for i in ids:
        num_referenced.append(track_references_antipodypedia(i))
    print('Successfully queried Antibodypedia')
    return links, num_referenced

#----------------- Compute Score -----------------#
def compute_score(row):
    '''
    Compute weighted score for protein selection based on article counts (favor fewer), 
    number of referenced antibodies, interactions GO score, and UV score.
    '''

    # late more important, mock less important ercc6 rated high by score, 
    # present at 3 timepoints, mock is one of them (fix uv score)
    # any term present gets increase score if found in abstract

    score = row['Total Article Count (normalized)']*-1 + row['Term Article Count (normalized)'] + row['Ref-ed Antibodies via Antibodypedia (normalized)'] + row['\'String\' Interactions (normalized)']*1.5 + row['GO Score (normalized)'] + row['UV Score (normalized)']*2
    return score

#----------------- Main -----------------#
def main(input_csv, output_csv):
    '''
    Read a list of proteins from a CSV file, query PubMed for the
    number of articles associated with each protein, string interactions, 
    query Antibodypedia.com for referenced antibodies, find associated 
    GO terms and compute a selection score, save results to CSV file.
    '''

    # Track start time
    start_time = time.time()

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
    links, num_referenced = query_antibodypedia(uniprots)
    total_article_counts, term_article_counts = query_pubmed(symbols)
    mapping = uniprot_to_string(uniprots)
    print('Not converted:', [u for u in uniprots if u not in mapping.keys()])
    strings = [mapping[uni] if uni in mapping.keys() else '' for uni in uniprots]
    interactions = query_string(strings)
    go_scores, go_terms = query_go_terms(uniprots)

    data_dict = {
        'Uniprot': uniprots,
        'Protein Symbol': symbols,
        'Total Article Count': total_article_counts,
        'Term Article Count': term_article_counts,
        '\'String\' Interactions': [interactions[string] for string in strings],
        'Antibodypedia Link': links,
        'Ref-ed Antibodies via Antibodypedia': num_referenced,
        'UV_treatment': uvs,
        'GO Score': go_scores,
        'GO Terms': go_terms
    }

    # Set up dataframe and normalize counts/scores
    df_all = pd.DataFrame(data_dict)
    df_all['UV Score'] = df_all['UV_treatment'].apply(lambda x: uv_score(x))
    for col in ['Ref-ed Antibodies via Antibodypedia', 'Total Article Count', 'Term Article Count', 
                '\'String\' Interactions', 'GO Score', 'UV Score']:
        df_all[col] =pd.to_numeric(df_all[col])
        df_all[f'{col} (normalized)'] = df_all[col]/df_all[col].std()

    # Compute weighted score for protein selection
    df_all['Overall Score'] = df_all.apply(lambda x: compute_score(x), axis=1)
    print('Computed scores')
    cols = ['Uniprot', 'Protein Symbol','Antibodypedia Link', 'Ref-ed Antibodies via Antibodypedia',
            'Total Article Count', 'Term Article Count', '\'String\' Interactions', 'GO Score', 
            'UV Score','Overall Score', 'GO Terms']
    df_all = df_all[cols]

    # Sort by overall score and save to CSV
    df_all.sort_values(by=['Overall Score'], ascending=False, inplace=True)
    df_all.to_csv(final_output_csv)

    # Track end time
    endtime = time.time()
    total_time = endtime - start_time
    print(f'Time taken: {total_time:.2f} seconds')
    print('Protein selection data saved to', final_output_csv)

# Main: Adjust input and output CSV file names accordingly on command line
if __name__ == "__main__":
    input_csv = sys.argv[1]
    final_output_csv = sys.argv[2]
    main(input_csv, final_output_csv)