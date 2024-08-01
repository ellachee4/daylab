# Author: Ella Chee
# Date: July 2024
# Description: Script to scrape antibody availability data from antibodypedia.com 
# for a list of UniProt IDs, given a CSV list of proteins, returns CSV containing
# the antibody link, number of antibodies, and number of providers for each protein.

# Imports
import sys
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

def scrape_antibodypedia_data(uniprot_id):
    '''Scrapes antibodypedia.com for the access link, number of antibodies and providers for a given UniProt ID'''

    # Set up the Selenium WebDriver
    driver = webdriver.Chrome()

    # Base URL
    base_url = 'https://www.antibodypedia.com/explore/uniprot%3A'

    # Add UniProt ID to base URL
    url = f'{base_url}{uniprot_id}'

    # Navigate to the URL
    driver.get(url)

    # Wait for the table to load
    WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.ID, "search_results_table"))
    )
    print('Found table for UNIPROT:', uniprot_id)

    data = []

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


    # Append the data to the list
    data.append({
        'UNPROT': uniprot_id,
        'Antibody Link': antibodies_link,
        'Number of Antibodies': number_of_antibodies,
        'Number of Providers': number_of_providers
    })

    driver.quit()
    return data

def main():
    '''Main function to scrape the data from antibodypedia.com and save it to a CSV file, 
        given a CSV of UniProt IDs'''

    # Read in CSV (change file name accordingly) to dataframe
    df = pd.read_csv('unique_to_mTbG4P_dedup.csv')
    uniprot_df = df[['UNIPROT']].copy()

    # Scrape data for each UniProt ID in dataframe
    for index, row in uniprot_df.iterrows():
        uniprot_id = row['UNIPROT']
        data = scrape_antibodypedia_data(uniprot_id)
        uniprot_df.loc[index, 'Antibody Link'] = data[0]['Antibody Link']
        uniprot_df.loc[index, 'Number of Antibodies'] = data[0]['Number of Antibodies']
        uniprot_df.loc[index, 'Number of Providers'] = data[0]['Number of Providers']
        print('Scraped data for UniProt ID:', uniprot_id)
    
    # Save the updated dataframe to a new CSV file (change file name as needed)
    uniprot_df.to_csv('antibody_availability_data.csv', index=False)

main()