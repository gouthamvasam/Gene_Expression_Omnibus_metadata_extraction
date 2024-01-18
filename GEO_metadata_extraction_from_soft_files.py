#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 11:10:11 2023

@author: gouthamvasam
"""
# GEO_metadata_extraction_from_soft_files

import os
import pandas as pd
import GEOparse
from Bio import Entrez
from Bio import Medline
import xmltodict
import requests
from bs4 import BeautifulSoup

os.chdir('path_to_working_directory')

# Read the GSE IDs from a text file
with open('desired_studies.txt', 'r') as file:
    gse_ids = [line.strip() for line in file]

metadata_path = 'path_to_save_soft_files'
output_file = 'path_to_save_output.tsv'

def get_series_metadata(gse_id, metadata_path):
    gse = GEOparse.get_GEO(geo=gse_id, destdir=metadata_path, silent=True)
    series_dfs = []
    series = pd.Series(gse.metadata).str.join(', ')
    series.index = 'series_' + series.index
    for gsm_name, gsm in gse.gsms.items():
        sample = pd.Series(gsm.metadata).str.join(', ')
        sample.index = 'sample_' + sample.index
        series_dfs.append(pd.concat((series, sample)))
    series_metadata = pd.concat(series_dfs, axis=1).T.set_index(['series_geo_accession', 'sample_geo_accession'])
    genome_build = series_metadata['sample_data_processing'].str.extract('Genome_build:(.+?),') 
    series_metadata['genome_build'] = genome_build
    return series_metadata.reset_index()

def get_paper_info(pubmed_id):
    handle = Entrez.efetch(db="pubmed", id=pubmed_id, retmode='xml', rettype='abstract')
    doc = xmltodict.parse(handle.read())
    try:
        doc = doc['PubmedArticleSet']['PubmedArticle']['MedlineCitation']['Article']
    except:
        pass
    paper_dict = {}
    try:
        paper_dict['year'] = doc['ArticleDate']['Year']
        paper_dict['month'] = doc['ArticleDate']['Month']
        paper_dict['day'] = doc['ArticleDate']['Day']
    except:
        pass
    try:
        paper_dict['title'] = doc['ArticleTitle']
    except:
        pass
    try:
        if isinstance(doc['Abstract']['AbstractText'], list) or isinstance(doc['Abstract']['AbstractText'], dict):
            abstract = [item['#text'] for item in doc['Abstract']['AbstractText'] if '#text' in item]
            paper_dict['abstract'] = ' '.join(abstract)
        else:
            paper_dict['abstract'] = doc['Abstract']['AbstractText']
    except:
        pass
    try:
        paper_dict['journal_title'] = doc['Journal']['Title']
    except:
        pass
    return paper_dict

def get_series_citation(gse_id):
    Entrez.email = "your.email@example.com"
    gse_url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}'
    gse_page_content = requests.get(gse_url).content
    soup = BeautifulSoup(gse_page_content, "html.parser")
    # Find the PMID
    citation_text_element = soup.find('td', text='Citation(s)')
    
    if citation_text_element is not None:
        citation_text_element = citation_text_element.find_next_sibling('td')
        citation_text = citation_text_element.get_text(strip=True, separator=' ')
        pmid = citation_text.split(' ')[-1]
        # Retrieve the citation information from PubMed
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
        record = Medline.read(handle)
        
        # Format the citation information
        citation = "{authors}. {title}. {journal} {date}. PMID: {pmid}".format(
            authors=", ".join(record["AU"]),
            title=record["TI"],
            journal=record["JT"],
            date=record["DP"],
            pmid=pmid
        )
        return citation
    else:
        return ''
    
# Create an empty DataFrame to store the metadata for all GSE IDs
all_metadata_df = pd.DataFrame()

# Loop over all GSE IDs and retrieve their metadata
for gse_id in gse_ids:
    # Get metadata for the current GSE ID
    metadata_df = get_series_metadata(gse_id, metadata_path)
    
    # Get citation for the current GSE ID
    citation = get_series_citation(gse_id)
    metadata_df['series_citation'] = citation

    # Check if the `series_pubmed_id` column is present in the metadata DataFrame
    if 'series_pubmed_id' in metadata_df.columns:
        # Add publication information to the metadata DataFrame
        pubmed_id = metadata_df.loc[0, 'series_pubmed_id']
        paper_info = get_paper_info(pubmed_id)
        for key, value in paper_info.items():
            if len(value) != len(metadata_df.index):
                value = pd.Series(value)
            metadata_df[key] = value

    # Append the metadata DataFrame to the all_metadata_df DataFrame
    all_metadata_df = pd.concat([all_metadata_df, metadata_df])

# Write the all_metadata_df DataFrame to a TSV file
all_metadata_df.to_csv(output_file, sep='\t', index=False)