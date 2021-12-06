#! /usr/bin/env python

import pandas as pd
import sys
import argparse
import re
import os
from datetime import date
# import subprocess
import glob
import shutil



def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("-p", '--path',  default = '', help="REQUIRED: path to metadata sheets")
    parser.add_argument("-s", '--submitter_name', default = '',  help = 'REQUIRED: submitter name for GISAID')

    options = parser.parse_args(args)
    return options

if __name__ == '__main__':
    
    # get user input and check that none are missing
    options = getOptions()  
    missing_flags = 0
    
    metadata_dir = options.path
    if metadata_dir == '':
        print('  ERROR must specify "-p or --path"')
        missing_flags =+ 1
    else:
        # set as working directory
        os.chdir(metadata_dir)

        
    submitter_name = options.submitter_name
    if submitter_name == '':
        print('  ERROR must specify "--submitter_name"')
        missing_flags =+ 1
    
    today_date = str(date.today())
    
    if missing_flags > 0:
        sys.exit()

    
    # print input variables:
    print('  INPUT FILES AND FLAGS:')
    print('')
    print('  metadata path: \n    %s\n' % metadata_dir)
    print('  submitter name: \n    %s\n' % submitter_name)
    print('  date:\n    %s\n' % today_date)
    print('')
    print('')
    
    
    # read in files
    ### from ncbi metadata-*-processed-ok.tsv
    for file in glob.glob('metadata-*-processed-ok.tsv'):                     
        ncbi_in = pd.read_csv(file, sep = '\t')
    #### pull out accession_id for joining
    for row in range(ncbi_in.shape[0]):
        filename = ncbi_in.filename[row]
        accession_id = filename.split('.')[0]
        if re.search('_', accession_id):
            accession_id = filename.split('_')[0]
        if re.search('CO-CDPHE-', accession_id):
            accession_id = accession_id.split('CO-CDPHE-')[1]
        ncbi_in.at[row, 'accession_id'] = accession_id
    ncbi_in = ncbi_in.set_index('accession_id')
    
    ### from gisaid 
    for file in glob.glob('gisaid_hcov-19_*.tsv'):
        gisaid_in = pd.read_csv(file, sep = '\t')
        gisaid_in = gisaid_in.rename(columns = {'Virus name' : 'virus_name'})
    ### pull out accession_id for joining
    for row in range(gisaid_in.shape[0]):
        virus_name = gisaid_in.virus_name[row]
        accession_id = re.findall('CO-CDPHE-([0-9\-\_]+)/', virus_name)[0]
        gisaid_in.at[row, 'accession_id'] = accession_id
    gisaid_in = gisaid_in.set_index('accession_id')
        
                          
    ncbi_gisaid_merged = ncbi_in.join(gisaid_in, how = 'left') 
#     ncbi_gisaid_merged = ncbi_gisaid_merged.reset_index()
    
    ### from genbank    
    ###### rename file genbank sends if not seqids.txt
    if os.path.exists('accessions.txt'):
        shutil.move('accessions.txt', 'seqids.txt')
    else:
        print('  GenBank file has expected file name\n')
        
    genbank = pd.read_csv('seqids.txt', sep = '\t', header = None, names = ['col1', 'GenBank'])
    for row in range(genbank.shape[0]):
        accession_id = re.findall('CO-CDPHE-([0-9\-\_]+)/', genbank.col1[row])[0]
        genbank.at[row, 'accession_id'] = accession_id
    genbank = genbank.set_index('accession_id')
    ncbi_gisaid_genbank_merged = ncbi_gisaid_merged.join(genbank, how = 'outer')
    
    ### fiiltered results file to get the seq_run number
    metadata_readin = pd.read_csv('filtered_results_subset_metadata.tsv', dtype = {'accession_id' : object}, sep = '\t')
    col_order = ['accession_id', 'seq_run']
    metadata_readin = metadata_readin[col_order]
    metadata_readin = metadata_readin.set_index('accession_id')


    all_merged = metadata_readin.join(ncbi_gisaid_genbank_merged, how = 'outer')
    all_merged = all_merged.reset_index()
    
    
    
    
    # filter for complete tab
    complete = all_merged[~all_merged.GenBank.isna()]
    complete = complete.reset_index(drop = True)
    
    ### add columns
    complete['submitter'] = submitter_name
    complete['Isolation Source'] = 'patient isolate'
    complete['Submission Date'] = today_date
    
    for row in range(complete.shape[0]):
        accession_id = complete.accession_id[row]
        isolate_name = 'CO-CDPHE-%s' % accession_id
        complete.at[row, 'isolate/sample_name'] = isolate_name

    
    rename_col = { 'virus_name':'Virus name', 
                  'title' :'sample_title', 
#                   'platform':'instrument_model', 
                  'bioproject_accession' :'BioProject', 
                  'biosample_accession': 'BioSample', 
                  'accession' :'SRA', 
                  'Accession ID' : 'GISAID'}
    
    completed_tab_columns = ['submitter', 'accession_id', 'seq_run', 'Virus name', 'isolate/sample_name',
       'sample_title', 'instrument_model', 'Isolation Source',
       'Collection date', 'Lineage', 'Clade', 'BioProject', 'BioSample', 'SRA',
       'GenBank', 'GISAID', 'Submission Date']
    complete = complete.rename(columns = rename_col)
    complete = complete[completed_tab_columns]
    complete = complete.sort_values(by = 'seq_run')
    
    outfile = 'COMPLETED_%s_metadata.tsv' % today_date
    complete.to_csv(outfile, sep = '\t', index = False)
    
    
    # filte for ongoing tab
    ongoing = all_merged[all_merged.GenBank.isna()]
    if ongoing.shape[0] > 0:
        ongoing = ongoing.reset_index(drop = True)

        ### add columns
        ongoing['submitter'] = submitter_name
        ongoing['Isolation Source'] = 'patient isolate'
        ongoing['Submission Date (genbank - final)'] = ''

        for row in range(ongoing.shape[0]):
            accession_id = ongoing.accession_id[row]
            isolate_name = 'CO-CDPHE-%s' % accession_id
            ongoing.at[row, 'isolate/sample_name'] = isolate_name


        rename_col = { 'virus_name':'Virus name', 
                      'title' :'sample_title', 
#                       'platform':'instrument_model', 
                      'bioproject_accession' :'BioProject', 
                      'biosample_accession': 'BioSample', 
                      'accession' :'SRA', 
                      'Accession ID' : 'GISAID',
                     'seq_run' : 'seq_run'}

        ongoing_tab_columns = ['submitter', 'accession_id', 'seq_run', 'Virus name',
           'isolate/sample_name', 'sample_title', 'instrument_model',
           'Isolation Source', 'Collection date', 'Lineage', 'Clade',
           'BioProject', 'BioSample', 'SRA', 'GenBank', 'GISAID',
           'Submission Date (genbank - final)']

        ongoing = ongoing.rename(columns = rename_col)
        ongoing = ongoing[ongoing_tab_columns]
        ongoing = ongoing.sort_values(by = 'seq_run')

        outfile = 'ONGOING_genbank_missing_%s_metadata.tsv' % today_date
        ongoing.to_csv(outfile, sep = '\t', index = False)
        
        print('%d sequences rejected by Genbank. Check %s file for details.' %(ongoing.shape[0], outfile))
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    