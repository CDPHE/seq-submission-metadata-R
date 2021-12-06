#! /usr/bin/env python

import pandas as pd
import sys
import argparse
import re
import os
from datetime import date
import subprocess
import glob

def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("-m", '--metadata',  default = '', help="REQUIRED: path to metadata sheets")

    options = parser.parse_args(args)
    return options
    
    
if __name__ == '__main__':
    
    # get user input and check that none are missing
    options = getOptions()  
    missing_flags = 0
    
    metadata_dir = options.metadata
    if metadata_dir == '':
        print('  ERROR must specify "--metadata"')
        missing_flags =+ 1
    else:
        # set as working directory
        os.chdir(metadata_dir)
   
    today_date = str(date.today())
    
    if missing_flags > 0:
        sys.exit()

    # get teh gisaid file input
    for file in glob.glob('gisaid_hcov-19_*.tsv'):
        gisaid_input_path = file
        
    # get the filtered data file input from the gisaid metadata formatting script
    for file in glob.glob('filtered_results_subset_metadata.tsv'):
        sample_df_path = file
    
    
    # print input variables:
    print('  *** INPUT FILES AND FLAGS:')
    print('')
    print('  metadata path: \n    %s\n' % metadata_dir)
    print('  gisaid processed path: \n    %s\n' % gisaid_input_path)
    print('  filtered sample list and metadata path: \n    %s\n' % sample_df_path)
    print('  date:\n    %s\n' % today_date)
    print('')
    print('')
    
        
    ##################
    ### read in the gisaid input data
    gisaid_processed_df = pd.read_csv(gisaid_input_path, sep = '\t') 
    
    rename_col = {'Virus name' : 'gisaid_virus_name',
                  'Accession ID' : 'gisaid_accession'}
    gisaid_processed_df = gisaid_processed_df.rename(columns = rename_col)
    
    ### get the co-cdphe-accession id from the covv_virus_anme column to use for merging
    def get_cdphe_accession(virus_name):
        accession = virus_name.split('/')[2]
        return accession
    
    gisaid_processed_df['cdphe_accession_id'] = gisaid_processed_df.apply(lambda x:get_cdphe_accession(x.gisaid_virus_name), axis = 1)
    
    ### filter the data set to only the columns we want
    keep_columns = ['cdphe_accession_id', 'gisaid_accession', 'gisaid_virus_name']
    gisaid_processed_df = gisaid_processed_df[keep_columns]
    
    ### set index for merging...
    gisaid_processed_df = gisaid_processed_df.set_index('cdphe_accession_id')
    
    
    #######################
    ### read in filtered sample dataset and then merge with gisaid data
    df = pd.read_csv(sample_df_path, sep = '\t', dtype = {'accession_id' : object})
    df = df.rename(columns = {'fasta_header' : 'cdphe_accession_id'})
    
    starting_col = ['cdphe_accession_id', 'collection_date']
    
    df = df.set_index('cdphe_accession_id')
    
    biosample_df = df.join(gisaid_processed_df, how = 'left')
    biosample_df = biosample_df.reset_index()
    
    
    
    ########################
    # create ncbi biosample file
    ### begin with starting with the joined dataset and renaming some columns
    col_rename = {'cdphe_accession_id' : 'sample_name'}
    biosample_df = biosample_df.rename(columns = col_rename)
    
    biosample_df['sample_title'] = 'PCR tiled amplicon WGS of SARS-CoV-2'
    biosample_df['bioproject_accession'] = 'PRJNA686984'
    biosample_df['organism'] = 'Severe acute respiratory syndrome coronavirus 2'
    biosample_df['collected_by'] = 'Colorado Department of Public Health and Environment'
    biosample_df['geo_loc_name'] = 'USA: Colorado'
    biosample_df['host'] = 'Homo sapiens'
    biosample_df['host_disease'] = 'severe accute respiratory syndrome'
    biosample_df['isolate'] = biosample_df.sample_name.tolist()
    biosample_df['isolation_source'] = 'patient isolate'
    biosample_df['sequenced_by'] = 'Colorado Department of Public Health and Environment'
    
    remaining_cols = ['antiviral_treatment_agent', 'collection_device', 'collection_method',
       'date_of_prior_antiviral_treat', 'date_of_prior_sars_cov_2_infection',
       'date_of_sars_cov_2_vaccination', 'exposure_event', 'geo_loc_exposure',
       'host_age','host_anatomical_material', 'host_anatomical_part', 'host_body_product',
       'host_disease_outcome', 'host_health_state', 'host_recent_travel_loc',
       'host_recent_travel_return_date', 'host_sex', 'host_specimen_voucher',
       'host_subject_id', 'lat_lon', 'passage_method', 'passage_number',
       'prior_sars_cov_2_antiviral_treat', 'prior_sars_cov_2_infection',
       'prior_sars_cov_2_vaccination', 'purpose_of_sampling',
       'purpose_of_sequencing', 'sars_cov_2_diag_gene_name_1',
       'sars_cov_2_diag_gene_name_2', 'sars_cov_2_diag_pcr_ct_value_1',
       'sars_cov_2_diag_pcr_ct_value_2', 'vaccine_received',
       'virus_isolate_of_prior_infection', 'description']
    
    for col in remaining_cols:
        biosample_df[col] = ''
        
    biosample_cols = ['sample_name', 'sample_title', 'bioproject_accession', 'organism',
       'collected_by', 'collection_date', 'geo_loc_name', 'host',
       'host_disease', 'isolate', 'isolation_source',
       'antiviral_treatment_agent', 'collection_device', 'collection_method',
       'date_of_prior_antiviral_treat', 'date_of_prior_sars_cov_2_infection',
       'date_of_sars_cov_2_vaccination', 'exposure_event', 'geo_loc_exposure',
       'gisaid_accession', 'gisaid_virus_name', 'host_age',
       'host_anatomical_material', 'host_anatomical_part', 'host_body_product',
       'host_disease_outcome', 'host_health_state', 'host_recent_travel_loc',
       'host_recent_travel_return_date', 'host_sex', 'host_specimen_voucher',
       'host_subject_id', 'lat_lon', 'passage_method', 'passage_number',
       'prior_sars_cov_2_antiviral_treat', 'prior_sars_cov_2_infection',
       'prior_sars_cov_2_vaccination', 'purpose_of_sampling',
       'purpose_of_sequencing', 'sars_cov_2_diag_gene_name_1',
       'sars_cov_2_diag_gene_name_2', 'sars_cov_2_diag_pcr_ct_value_1',
       'sars_cov_2_diag_pcr_ct_value_2', 'sequenced_by', 'vaccine_received',
       'virus_isolate_of_prior_infection', 'description']
    

    biosample_df = biosample_df[biosample_cols]
    
    
    #### check that all samples have gisaid_accession...
    crit = biosample_df.gisaid_accession.isna()
    missing_df = biosample_df[crit]
    if missing_df.shape[0] > 0:
        print('\n  WARNING: the following sequences are missing GISAID information.')
        print('  Check that gisaid submissions have been processed')
        print(missing_df.sample_name)
        print('')
    
    # save biosample df
    outfile = os.path.join(metadata_dir, 'ncbi_biosample_submission_%s_metadata.tsv' % today_date )
    biosample_df.to_csv(outfile, index = False, sep = '\t')
    print('  *** WRITING BIOSAMPLE METADATA FILE TO:')
    print('    %s\n' % outfile)
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    