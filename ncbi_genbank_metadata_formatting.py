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
    parser.add_argument("-p", '--path',  default = '', help="REQUIRED: path to metadata sheets")
    parser.add_argument("-t", '--seq_type', default = '',  help = "REQUIRED: sequencing type as a string in quotes that has to match 'Illumina MiSeq', 'Illumina NextSeq', 'Illumina NovaSeq', 'Oxford Nanopore GridION'")

    options = parser.parse_args(args)
    return options

if __name__ == '__main__':
    
    # get user input and check that none are missing
    options = getOptions()  
    missing_flags = 0
    
    metadata_dir = options.path
    if metadata_dir == '':
        print('  ERROR must specify "--metadata"')
        missing_flags =+ 1
    else:
        # set as working directory
        os.chdir(metadata_dir)
        
    tech_platform = options.seq_type
    if tech_platform == '':
        print('  ERROR must specify "--seq_type"')
        missing_flags =+ 1
    
    today_date = str(date.today())
    
    if missing_flags > 0:
        sys.exit()

    
    ## get the gisaid and sra input file paths
    for file in glob.glob('metadata-*-processed-ok.tsv'):
        sra_input_path = file
    
    for file in glob.glob('gisaid_hcov-19_*.tsv'):
        gisaid_input_path = file
        
    # print input variables:
    print('  *** INPUT FILES AND FLAGS:')
    print('')
    print('  metadata path: \n    %s\n' % metadata_dir)
    print('  sequencing platform: \n    %s\n' % tech_platform)
    print('  date:\n    %s\n' % today_date)
    print('  metadata-processed-ok.tsv path: \n    %s\n' % sra_input_path)
    print('  gisaid processed path: \n    %s\n' % gisaid_input_path)
    print('')
    print('')
    
    # error out if tech platform does not match a given option
    tech_platform_options = ['Illumina MiSeq', 'Illumina NextSeq', 'Oxford Nanopore GridION', 'Illumina NovaSeq']
    if tech_platform not in tech_platform_options:
        print('')
        print('  ERROR: invalid input for "--seq_type"')
        print('    choose among one of the following options:')
        for platform in tech_platform_options:
            print('      %s' % platform)
        print('')
        sys.exit()
    
    
    ##################
    ### read in the gisaid and sra input data
    gisaid_processed_df = pd.read_csv(gisaid_input_path, sep = '\t')
    sra_processed_df = pd.read_csv(sra_input_path, sep = '\t')
    
    
    ###################
    ### create the ncbi genbank dataframe
    ### begin with the sra_processed_df
    def get_cdphe_accession_id_covmin(filename):
        cdphe_accession = filename.split('.')[0]
        return cdphe_accession
    
    def get_cdphe_accession_id_nexseq(filename):
        accession_id = filename.split('.')[0]
        cdphe_accession = 'CO-CDPHE-%s' % accession_id
        return cdphe_accession
    
    def get_cdphe_accession_id_miseq(filename):
        accession_id = filename.split('_')[0]
        cdphe_accession = 'CO-CDPHE-%s' % accession_id
        return cdphe_accession
        
    starting_columns = ['filename', 'biosample_accession', 'accession', ]
    rename_columns = {'biosample_accession' : 'BioSample',
                      'accession' : 'SRA'}
    
    ncbi_genbank_df = sra_processed_df[starting_columns]
    ncbi_genbank_df = ncbi_genbank_df.rename(columns = rename_columns)
    
    ## craete cdphe accession column for merging
    if tech_platform == 'Oxford Nanopore GridION':
        ncbi_genbank_df['cdphe_accession'] = ncbi_genbank_df.apply(lambda x:get_cdphe_accession_id_covmin(x.filename), axis = 1)
    elif tech_platform == 'Illumina NextSeq':
        ncbi_genbank_df['cdphe_accession'] = ncbi_genbank_df.apply(lambda x:get_cdphe_accession_id_nexseq(x.filename), axis = 1)
    elif tech_platform == 'Illumina MiSeq' :
        ncbi_genbank_df['cdphe_accession'] = ncbi_genbank_df.apply(lambda x:get_cdphe_accession_id_miseq(x.filename), axis = 1)
#     elif tech_platform == 'Illumina NovaSeq':
#         ncbi_genbank_df['cdphe_accession'] = ncbi_genbank_df.apply(lambda x:get_cdphe_accession_id_miseq(x.filename), axis = 1)
    
    ### merge gisaid columns
    ### prep gisaid columns
    def get_accession_from_gisaid_virus_name(virus_name):
        cdphe_accession = virus_name.split('/')[2]
        return cdphe_accession
        
    keep_columns = ['Virus name', 'Accession ID', 'Collection date']
    rename_col = {'Virus name' : 'Sequence_ID',
                  'Collection date' : 'collection-date'}
    gisaid_processed_df = gisaid_processed_df[keep_columns]
    gisaid_processed_df = gisaid_processed_df.rename(columns = rename_col)
    
    gisaid_processed_df['cdphe_accession'] = gisaid_processed_df.apply(lambda x:get_accession_from_gisaid_virus_name(x.Sequence_ID), axis = 1)
    
    ### merge
    ncbi_genbank_df = ncbi_genbank_df.set_index('cdphe_accession')
    gisaid_processed_df = gisaid_processed_df.set_index('cdphe_accession')
    
    ncbi_genbank_df = ncbi_genbank_df.join(gisaid_processed_df, how = 'left')
    ncbi_genbank_df = ncbi_genbank_df.reset_index()
    
    ## rename some columns
    col_rename = {'Accession ID' : 'gisaid_id',
                  'cdphe_accession' : 'isolate'}
    ncbi_genbank_df = ncbi_genbank_df.rename(columns = col_rename)

    
    ### create remaining columns
    def create_ncbi_genbank_note(gisaid_id):
        if gisaid_id:
            note = 'GISAID_ID: %s' % gisaid_id
            return note
        else:
            return ''
        
    ncbi_genbank_df['country'] = 'USA'
    ncbi_genbank_df['host'] = 'Homo sapiens'
    ncbi_genbank_df['isolation-source'] = 'patient isolate'
    ncbi_genbank_df['note'] = ncbi_genbank_df.apply(lambda x:create_ncbi_genbank_note(x.gisaid_id), axis = 1)
                       
    ncbi_genbank_columns = ['Sequence_ID', 'country', 'host', 'isolate', 'collection-date', 'isolation-source', 
                            'BioSample', 'SRA', 'note']

    ncbi_genbank_df = ncbi_genbank_df[ncbi_genbank_columns]
    ncbi_genbank_df = ncbi_genbank_df.reset_index(drop = True)
    
    ### check for missing gisaid ids
    crit = ncbi_genbank_df.Sequence_ID.isna()
    missing_df = ncbi_genbank_df[crit]
    if missing_df.shape[0] > 0:
        print('\n  WARNING: the following sequences are missing GISAID information.')
        print(missing_df.isolate)
        print('')

    ### save ncbi_genbank_df
    outfile = os.path.join(metadata_dir, 'ncbi_genbank_submission_%s_metadata.tsv' % today_date)
    ncbi_genbank_df.to_csv(outfile, sep = '\t', index = False)
    print('\n  WRITING ncbi genbank metdata to:')
    print('    %s' % outfile)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    