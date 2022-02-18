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
    
    print('')
    print('  ***** RUNNING ncbi_sra_metadata_formatting.py *****')
    print('  last updated: 2022-02-18')
    print('  updates: fixed NextSeq 550 instrument model')
    print('')
    print('')
    
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

    ## get the biosample data path
    for file in glob.glob('BioSampleObjects.txt'):
        biosample_input_path = file  
        
    ## get the fiiltered samples
    for file in glob.glob('filtered_results_subset_metadata.tsv'):
        filtered_data_path = file
        
        
    # print input variables:
    print('  *** INPUT FILES AND FLAGS:')
    print('')
    print('  metadata path: \n    %s\n' % metadata_dir)
    print('  sequencing platform: \n    %s\n' % tech_platform)
    print('  biosample_input_path: \n    %s\n' % biosample_input_path)
    print('  filtered_metadata_path: \n    %s\n' % filtered_data_path)
    print('  date:\n    %s\n' % today_date)
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

    ######## read in teh filtered data (we will need to pull the primer info)
    df = pd.read_csv(filtered_data_path, sep = '\t', dtype = {'accession_id' : object})
    col_rename = {'fasta_header' : 'cdphe_accession_id'}
    df = df.rename(columns = col_rename)
    
    col_keep = ['cdphe_accession_id', 'primer_set']
    df = df[col_keep]
    
    df = df.set_index('cdphe_accession_id')
        
    ########## read in teh biosample data
    biosample = pd.read_csv(biosample_input_path, sep = '\t')
    col_rename = {'BioProject' : 'bioproject_accession', 
                  'Accession' : 'biosample_accession',
                 'Sample Name' : 'cdphe_accession_id'}
    biosample = biosample.rename(columns = col_rename)
    col_keep = ['bioproject_accession', 'biosample_accession', 'cdphe_accession_id']
    biosample = biosample[col_keep]
    
    biosample = biosample.set_index('cdphe_accession_id')
    
    ###### join the df and the biosample to start the sra datatable
    sra = df.join(biosample, how = 'outer')
    sra = sra.reset_index()
    num_samples = sra.shape[0]
    
    if sra.shape[0] != df.shape[0] and sra.shape[0] != biosample.shape[0]:
        print('  WARNING: the biosample and sra datasets are not the same size')
    
    ## add in the additional columns
    def get_design_description(primer_set):
        design_description = 'Whole genome sequencing of SARS-CoV-2 using the %s protocol' % primer_set
        return design_description
    
    def get_fastq_name_ONT(cdphe_accession_id):
        fastq_name = '%s.fastq.gz' % cdphe_accession_id
        return fastq_name
        
    def get_fastq_name_ILLUMINA_single(cdphe_accession_id):
        accession_id = cdphe_accession_id.split('CO-CDPHE-')[1]
        fastq_name = '%s.fastq.gz' % accession_id
        return fastq_name
        
    def get_fastq_name_ILLUMINA_paired_R1(cdphe_accession_id):
        accession_id = cdphe_accession_id.split('CO-CDPHE-')[1]
        fastq_R1_name = '%s_R1_001.fastq.gz' % accession_id
        return fastq_R1_name
    
    def get_fastq_name_ILLUMINA_paired_R2(cdphe_accession_id):
        accession_id = cdphe_accession_id.split('CO-CDPHE-')[1]
        fastq_R2_name = '%s_R2_001.fastq.gz' % accession_id
        return fastq_R2_name
    
    sra['library_ID'] = range(num_samples)
    sra['title'] = 'PCR tiled amplicon WGS of SARS-CoV-2'
    sra['library_strategy'] = 'AMPLICON'
    sra['library_source'] = 'VIRAL RNA'
    sra['library_selection'] = 'PCR'
    sra['design_description'] = sra.apply(lambda x:get_design_description(x.primer_set), axis = 1)
    sra['filetype'] = 'fastq'
    
    if tech_platform in ['Illumina NextSeq']:
        sra['library_layout'] = 'single'
        sra['platform'] = 'ILLUMINA'
        sra['instrument_model'] = 'NextSeq 550'
        sra['filename'] = sra.apply(lambda x:get_fastq_name_ILLUMINA_single(x.cdphe_accession_id), axis = 1)
        sra['filename2'] = ''
        
    elif tech_platform in ['Illumina MiSeq']:
        sra['library_layout'] = 'paired'
        sra['platform'] = 'ILLUMINA'
        sra['instrument_model'] = 'Illumina MiSeq'
        sra['filename'] = sra.apply(lambda x:get_fastq_name_ILLUMINA_paired_R1(x.cdphe_accession_id), axis = 1)
        sra['filename2'] = sra.apply(lambda x:get_fastq_name_ILLUMINA_paired_R2(x.cdphe_accession_id), axis = 1)
        

    elif tech_platform in ['Illumina NovaSeq']:
        sra['library_layout'] = 'single'
        sra['platform'] = 'ILLUMINA'
        sra['instrument_model'] = 'Illumina NovaSeq 6000'
        sra['filename'] = sra.apply(lambda x:get_fastq_name_ILLUMINA_single(x.cdphe_accession_id), axis = 1)
        sra['filename2'] = ''
        
    elif tech_platform in ['Oxford Nanopore GridION']:
        sra['library_layout'] = 'single'
        sra['platform'] = 'OXFORD_NANOPORE'
        sra['instrument_model'] = 'GridION'
        sra['filename'] = sra.apply(lambda x:get_fastq_name_ONT(x.cdphe_accession_id), axis = 1)
        sra['filename2'] = ''
        
    #### remove the unneeded columns and get columns in correct order      
    col_order = ['bioproject_accession', 'biosample_accession', 'library_ID', 'title',
       'library_strategy', 'library_source', 'library_selection',
       'library_layout', 'platform', 'instrument_model', 'design_description',
       'filetype', 'filename', 'filename2']
    
    sra = sra[col_order]
    
    ### save sra df
    outfile = os.path.join(metadata_dir, 'ncbi_sra_submission_%s_metadata.tsv' % today_date)
    sra.to_csv(outfile, sep = '\t', index = False)
    print('\n  WRITING ncbi sra metdata to:')
    print('    %s' % outfile)
    

    
