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
    parser.add_argument("-f", "--fasta_name", default = '', help="REQUIRED: fasta file name to upload to GISAID")
    parser.add_argument("-s", '--submitter_name', default = '',  help = 'REQUIRED: submitter name for GISAID')
    parser.add_argument("-r", '--rerun_file', default = '',  help = 'REQUIRED: path and file name for text file with rerun information as downloaded from Google Drive')
    parser.add_argument("-c", '--completed_file', default = '',  help = 'REQUIRED: path and file name for text file with completed submission information as downloaded from Google Drive')
    parser.add_argument("-o", '--ongoing_file', default = '',  help = 'REQUIRED: path and file name for text file with ongoing submission information as downloaded from Google Drive')
    parser.add_argument('-z', '--coverage_min', default = 50,  type = int, help = 'OPTIONAL: minmium coverage')
    parser.add_argument("-y", '--coverage_max', default = 100, type = int,  help = 'OPTIONAL: maximum coverage')
    parser.add_argument("-t", '--seq_type', default = '',  help = "REQUIRED: sequencing type as a string in quotes that has to match 'Illumina MiSeq', 'Illumina NextSeq', 'Illumina NovaSeq', 'Oxford Nanopore GridION'")
    parser.add_argument('-p', '--primer_set', default = 'Artic V3')

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
        
    gisaid_fasta_name = options.fasta_name
    if gisaid_fasta_name == '':
        print('  ERROR must specify "--fasta_name"')
        missing_flags =+ 1
        
    submitter_name = options.submitter_name
    if submitter_name == '':
        print('  ERROR must specify "--submitter_name"')
        missing_flags =+ 1
        
    tech_platform = options.seq_type
    if tech_platform == '':
        print('  ERROR must specify "--seq_type"')
        missing_flags =+ 1
        
    rerun_file_path = options.rerun_file
    if rerun_file_path == '':
        print('  ERROR must specify "--rerun_file"')
        missing_flags =+ 1
        
    completed_submissions_path = options.completed_file
    if completed_submissions_path == '':
        print('  ERROR must specify "--completed_file"')
        missing_flags =+ 1
        
    ongoing_submissions_path = options.ongoing_file
    if ongoing_submissions_path == '':
        print('  ERROR must specify "--ongoing_file"')
        missing_flags =+ 1
        
    primer_set = options.primer_set
    coverage_min = options.coverage_min
    coverage_max = options.coverage_max
    
    today_date = str(date.today())
    
    if missing_flags > 0:
        sys.exit()

    
    # print input variables:
    print('  *** INPUT FILES AND FLAGS:')
    print('')
    print('  metadata path: \n    %s\n' % metadata_dir)
    print('  name of fasta file to submit to GISAID: \n    %s\n' % gisaid_fasta_name)
    print('  submitter name for GISAID: \n    %s\n' % submitter_name)
    print('  sequencing platform: \n    %s\n' % tech_platform)
    print('  name of rerun file: \n    %s\n' % rerun_file_path)
    print('  name of completed submissions file: \n    %s\n' % completed_submissions_path)
    print('  name of ongoing submissions file: \n    %s\n' % ongoing_submissions_path)
    print('  min coverage cutoff: \n    %d\n' % coverage_min)
    print('  max coverage cutoff:\n    %d\n' % coverage_max)
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

    # Reading in all metadata sheets, subsetting out the columns needed, getting rid of blanks and controls, and pulling samples with >= 50% coverage
    df_list = []
    for file in glob.glob('results_*.tsv'):
        df = pd.read_csv(file, sep = '\t', dtype = {'accession_id' : object}, na_values = ['', 'sample failed assembly'])
        df_list.append(df)
    metadata_readin = pd.concat(df_list)
    metadata_readin = metadata_readin.reset_index(drop = True)
    metadata_readin = metadata_readin.rename(columns = {'percent_non_ambigous_bases' : 'coverage'})
    
    # check for surveillance data
    if 'surveillance' in metadata_readin.columns:
        surveillance = True
        print('  *** SURVEILLANCE:')
        print('     results.tsv files have surveillance data')
        print('\n')
        
    else:
        surveillance = False
        metadata_readin['surveillance'] = ''
        print('  *** SURVEILLANCE:')
        print('    results.tsv files DO NOT have surveillance data')
        print('\n')
    
    # check for primer data
    if 'primer_set' in metadata_readin.columns:
        primer_set = True
        primers_used = metadata_readin.primer_set.unique().tolist()
        
        # if COVIDSEQ 3 primers used update to Artic V3 primers
        for row in range(metadata_readin.shape[0]):
            primer_set = metadata_readin.primer_set[row]
            if primer_set == 'COVIDSEQ V3' or primer_set == 'COVIDSeqV3' or primer_set == 'COVIDSeq V3' or primer_set == 'V3':
                metadata_readin.at[row, 'primer_set'] = 'Artic V3'
        # fix Artic V4 primers....
        for row in range(metadata_readin.shape[0]):
            primer_set = metadata_readin.primer_set[row]
            if primer_set == 'V4' :
                metadata_readin.at[row, 'primer_set'] = 'Artic V4'        
        # if missing fill in with 'Artic V3'
        crit = metadata_readin.primer_set.isna()
        missing_primers = metadata_readin[crit]
        if missing_primers.shape[0] > 0 :
            seq_runs_with_missing_primers = missing_primers.seq_run.unique().tolist()
            metadata_readin.primer_set = metadata_readin.primer_set.fillna('Artic V3')
        else:
            seq_runs_with_missing_primers = []
        
        ### now continue onward
        primers_used_dict = dict(zip(metadata_readin.seq_run, metadata_readin.primer_set))    
        print('  *** PRIMER_SET:')
        print('     results.tsv file include primer information')
        if 'COVIDSEQ V3' in primers_used:
            print('     COVIDSEQ V3 primers will be renamed to Artic V3')

        if len(seq_runs_with_missing_primers) > 0:
            print('     The following seq runs did not include primer info will default to "Artic V3"')
            for seq_run in seq_runs_with_missing_primers:
                print('          %s' % seq_run)
                
        print('')
        print('     primers used include:')        
        for seq_run in primers_used_dict:
            print('          %s : %s' % (seq_run, primers_used_dict[seq_run]))
        print('\n')
        
    else:
        primer_set = 'Artic V3'
        metadata_readin['primer_set'] = primer_set
        primers_used_dict = dict(zip(metadata_readin.seq_run, metadata_readin.primer_set))
        
        print('  *** PRIMER_SET:')
        print('     results.tsv file DO NOT INCLUDE primer information\n')
        print('     will use %s as default; if this is not correct, specify primer set using "--primer_set"' % primer_set)
        for seq_run in primers_used_dict:
            print('          %s : %s' % (seq_run, primers_used_dict[seq_run]))
        print('\n')
     
    
    
    ##################
    # filter based on percent coverage and drop pos/ neg controls and drop unused columns... 
    crit1 = metadata_readin.coverage >= coverage_min
    crit2 = metadata_readin.coverage < coverage_max
    crit3 = ~metadata_readin.accession_id.str.contains('Blank|NC|NTC|ExtractionPositive|POS')   

    if coverage_max == 100:
        metadata_readin = metadata_readin[crit1 & crit3]
    else:
        metadata_readin = metadata_readin[crit1 & crit2 & crit3]

    col_order = ['accession_id', 'coverage', 'collection_date', 'fasta_header', 'seq_run', 'surveillance', 'primer_set']
    metadata_readin = metadata_readin[col_order]
    metadata_readin = metadata_readin.reset_index(drop = True)
    metadata_readin = metadata_readin.sort_values(by = 'coverage', ascending = False)
#     metadata_readin.collection_date = pd.to_datetime(metadata_readin.collection_date)
#     print(metadata_readin.head())
    
    
    ###############################  
    # rerun samples
    ### read in rerun samples
    rerun_df = pd.read_csv(rerun_file_path, sep = '\t', dtype = {'accession_id' : object})
    
    ##### fix first run names....
    for row in range(rerun_df.shape[0]):
        first_run = rerun_df.first_run[row]
        if not pd.isnull(first_run):
            if re.search('\(', first_run):
                rerun_df.at[row, 'first_run'] = first_run.split('(')[0]

    rerun_accessions = rerun_df.accession_id.unique().tolist()
    metadata_readin_accessions = metadata_readin.accession_id.unique().tolist()
    metadata_readin_seq_runs = metadata_readin.seq_run.unique().tolist()
    
    
    
    #### print to screen the list of rerun samples and also save to outfile
    crit = rerun_df.accession_id.isin(metadata_readin_accessions)
    temp =rerun_df[crit]
    
    if temp.shape[0] > 0:
        rerun_first_runs = temp.first_run.unique().tolist()
        rerun_re_runs = temp.re_run.unique().tolist()
        
        outfile = os.path.join(metadata_dir, 'samples_rerun_%s.tsv' % today_date)
        temp.to_csv(outfile, sep = '\t', index= False)
        
        print('  *** BAD NOODLES:')
        print('  WARNING: The samples below have been submitted under multiple projects.')
        print('           If you are processing the project in "first_run" these samples will be removed from the dataset:\n')
        print(temp[['accession_id', 'first_run', 're_run']])
        print('\n')

        #### drop the accessions if the re-run seq_run doesn't match the seq run in the metadata_readin
        to_drop_accessions_dict = dict(zip(temp.accession_id, temp.first_run))
#         print(to_drop_accessions_dict)

        rows_to_drop = []
        for row in range(metadata_readin.shape[0]):
            accession_id = metadata_readin.accession_id[row]
            seq_run = metadata_readin.seq_run[row]
            if accession_id in to_drop_accessions_dict:
                if seq_run == to_drop_accessions_dict[accession_id]:
                    # drop that row otherwise keep
                    rows_to_drop.append(row)


        drop_crit = ~metadata_readin.index.isin(rows_to_drop)
        metadata_readin = metadata_readin[drop_crit]
        metadata_readin = metadata_readin.reset_index(drop = True)
    
    
    #######################
    # check for duplicate accession ids
    dups = []
    accession_ids = []
    for row in range(metadata_readin.shape[0]):
        accession_id = metadata_readin.accession_id[row]
        if accession_id in accession_ids:
            dups.append(accession_id)
        else:
            accession_ids.append(accession_id)
            
    if len(dups) > 0:
        print('  *** WITHIN SUBMISSION BATCH DUPS:')
        print('    WARNING: the follwing samples were duplicated in the batch. The sample with higher coverage will be retained.\n')
        for dup in dups:
            print('          %s' % dup)
        print('\n')
    
    ##### drop duplicates
    metadata_readin = metadata_readin.drop_duplicates(subset = 'accession_id', keep = 'first')
        
        
    ################################    
    # check for wrong collection dates
    metadata_readin.collection_date = pd.to_datetime(metadata_readin.collection_date)
    crit1 = metadata_readin.collection_date.isna()
    crit2 = metadata_readin.collection_date < '2020-01-01' #date(2020, 1, 1)
    temp = metadata_readin[crit1 | crit2]
    if temp.shape[0] > 0:
        outfile = os.path.join(metadata_dir, 'samples_with_wrong_collection_date_%s.tsv' % today_date)
        temp.to_csv(outfile, sep = '\t', index= False)
        
        print('  *** WRONG/MISSING COLLECTION DATE:')
        print('    WARNING: The samples below are either missing a collection date or have the wrong collection date')
        print('             Cannot continue with submission until these dates are fixed\n')
        print(temp[['accession_id', 'seq_run', 'collection_date']])
        print('\n')

        
    ###############################  
    # completed submission
    ### read in complted submissions
    completed_submissions = pd.read_csv(completed_submissions_path, sep = '\t', dtype = {'accession_id' : object})
    completed_accessions = completed_submissions.accession_id.unique().tolist()
    metadata_readin_accessions = metadata_readin.accession_id.unique().tolist()
    
    #### drop those accessions that have already be submitted
    crit1 = ~metadata_readin.accession_id.isin(completed_accessions)
    metadata_readin = metadata_readin[crit1]
    metadata_readin = metadata_readin.reset_index(drop = True)
    
    #### print to screen those that have already been submitted and save to file
    crit2 = completed_submissions.accession_id.isin(metadata_readin_accessions)
    temp = completed_submissions[crit2]
    outfile = os.path.join(metadata_dir, 'samples_already_submitted_%s.tsv' % today_date)
    temp.to_csv(outfile, sep = '\t', index = False)
    
    if temp.shape[0] > 0:
        print('  *** ALREADY SUBMITTED:')
        print('  WARNING: The samples below have already been submitted (completed submissions) and are being deleted from the dataset:\n')
        print(temp[['submitter', 'accession_id', 'seq_run', 'Submission Date']])
        print('\n')

        
        
    ###############################  
    # ongoing submission
    ### read in ongoing submissions
    ongoing_submissions = pd.read_csv(ongoing_submissions_path, sep = '\t', dtype = {'accession_id' : object})
    ongoing_accessions = ongoing_submissions.accession_id.unique().tolist()
    metadata_readin_accessions = metadata_readin.accession_id.unique().tolist()
    
    #### drop those accessions that have already be submitted
    crit1 = ~metadata_readin.accession_id.isin(ongoing_accessions)
    metadata_readin = metadata_readin[crit1]
    metadata_readin = metadata_readin.reset_index(drop = True)
    
    #### print to screen those that have already been submitted and save to file
    crit2 = ongoing_submissions.accession_id.isin(metadata_readin_accessions)
    temp = ongoing_submissions[crit2]
    outfile = os.path.join(metadata_dir, 'samples_in_progress_already_%s.tsv' % today_date)
    temp.to_csv(outfile, sep = '\t', index = False)
    
    if temp.shape[0] > 0:
        print('  *** ONGOING SUBMISSIONS:')
        print('  WARNING: The samples below are in progress of already beiing submitted (ongoing submissions)')
        print('           and will be removed from this submission:\n')
        print(temp[['submitter', 'accession_id', 'Seq_run']])
        print('\n')

    


    
    ###################
    # to screen the number of rows of data
    print('  *** DONE FILTERING METADATA')
    print('    there are %d unique accessions to submit in this batch\n' % len(metadata_readin.accession_id.unique().tolist()))
    
    ####################
    # save the filtered metadata_readin df
    outfile = os.path.join(metadata_dir, 'filtered_results_subset_metadata.tsv' )
    metadata_readin.to_csv(outfile, index = False, sep = '\t')
    print('  *** WRITING FILTERED METADATA FILE TO:')
    print('    %s\n' % outfile)
    
    
    ####################
    # create gisaid metadata dataframe
    ### first begin with metadata_readin in data and do some manipulation with those columns
    starting_col = ['accession_id', 'collection_date', 'fasta_header', 'surveillance', 'primer_set']
    gisaid_df = metadata_readin[starting_col]
    
    col_rename = {'collection_date' : 'covv_collection_date',
                 'surveillance' : 'covv_sampling_strategy'}
    gisaid_df = gisaid_df.rename(columns = col_rename)
    
    def create_covv_virus_name(covv_collection_date, fasta_header):
        collection_year = str(covv_collection_date).split('-')[0]
        covv_virus_name = 'hCoV-19/USA/%s/%s' % (fasta_header, collection_year)
        return covv_virus_name
       
    gisaid_df['submitter'] = submitter_name
    gisaid_df['fn'] = '%s.fasta' % gisaid_fasta_name
    gisaid_df['covv_virus_name'] = gisaid_df.apply(lambda x:create_covv_virus_name(x.covv_collection_date, x.fasta_header), axis=1)
    gisaid_df['covv_type'] = 'betacoronavirus'
    gisaid_df['covv_passage'] = 'Original'
    gisaid_df['covv_location'] = 'North America / USA / Colorado'
    gisaid_df['covv_add_location'] = ''
    gisaid_df['covv_host'] = 'Human'
    gisaid_df['covv_add_host_info'] = ''
    gisaid_df['covv_gender'] = 'unknown'
    gisaid_df['covv_patient_age'] = 'unknown'
    gisaid_df['covv_patient_status'] = 'unknown'
    gisaid_df['covv_specimen'] = ''
    gisaid_df['covv_outbreak'] = ''
    gisaid_df['covv_last_vaccinated'] = ''
    gisaid_df['covv_treatment'] = ''
    gisaid_df['covv_seq_technology'] = tech_platform
    gisaid_df['covv_assembly_method'] = ''
    gisaid_df['covv_coverage'] = ''
    gisaid_df['covv_orig_lab'] = 'Colorado Department of Public Health and Environment'
    gisaid_df['covv_orig_lab_addr'] = '8100 E. Lowry Blvd. Denver CO, 80230'
    gisaid_df['covv_provider'] = ''
    gisaid_df['covv_subm_lab'] = 'Colorado Department of Public Health and Environment'
    gisaid_df['covv_subm_lab_addr'] = '8100 E. Lowry Blvd. Denver CO, 80230'
    gisaid_df['covv_subm_sample_id'] = ''
    gisaid_df['covv_authors'] = 'Laura Bankers, Molly C. Hetherington-Rauth, Diana Ir, Alexandria Rossheim, Michael Martin, Mandy Waters, Arianna Smith, Shannon R. Matzinger, Emily A. Travanty'
    gisaid_df['covv_comment'] = ''
    gisaid_df['comment_type'] = ''
    gisaid_df['covv_add_host_info'] = ''
    gisaid_df['covv_provider_sample_id'] = ''
    gisaid_df['covv_specimen'] = ''
    
    gisaid_cols = ['submitter', 'fn', 'covv_virus_name', 'covv_type', 'covv_passage',
       'covv_collection_date', 'covv_location', 'covv_add_location',
       'covv_host', 'covv_add_host_info', 'covv_sampling_strategy',
       'covv_gender', 'covv_patient_age', 'covv_patient_status',
       'covv_specimen', 'covv_outbreak', 'covv_last_vaccinated',
       'covv_treatment', 'covv_seq_technology', 'covv_assembly_method',
       'covv_coverage', 'covv_orig_lab', 'covv_orig_lab_addr',
       'covv_provider_sample_id', 'covv_subm_lab', 'covv_subm_lab_addr',
       'covv_subm_sample_id', 'covv_authors', 'covv_comment', 'comment_type']

    gisaid_df = gisaid_df[gisaid_cols]
    
    #### write out gisaid_df
    outfile = os.path.join(metadata_dir, 'gisaid_submission_%s_metadata.csv' % today_date)
    gisaid_df.to_csv(outfile, index = False)
    print('  *** WRITING GISAID METADATA FILE TO:')
    print('    %s\n' % outfile)
    
    
    ##################
    # create rename fasta file...
    ##### begin with columns from metadata_readin df
    starting_col = ['fasta_header', 'seq_run', 'collection_date', 'surveillance']
    rename_fasta = metadata_readin[starting_col]
    
    col_rename = {'fasta_header' : 'accession_id',
                  'seq_run' : 'proj_num'}   
    rename_fasta = rename_fasta.rename(columns = col_rename)
    
    rename_fasta['gisaid_name'] = rename_fasta.apply(lambda x:create_covv_virus_name(x.collection_date, x.accession_id), axis=1)
    
    def create_ncbi_name(surveillance, gisaid_name):
        if surveillance == 'baseline surveillance':
            ncbi_name = '%s [PRJNA686984] [keyword=purposeofsampling:baselinesurveillance]' % gisaid_name
        elif surveillance == 'targeted surveillance':
            ncbi_name = '%s [PRJNA686984] [keyword=purposeofsampling:targetedsurveillance]' % gisaid_name
        else:
            ncbi_name = '%s [PRJNA686984]' % gisaid_name
        return ncbi_name
        
    rename_fasta['ncbi_name'] = rename_fasta.apply(lambda x:create_ncbi_name(x.surveillance, x.gisaid_name), axis = 1)
    
    col_order = ['accession_id', 'gisaid_name', 'proj_num', 'ncbi_name']
    rename_fasta = rename_fasta[col_order]
    
    #### save to file
    outfile = os.path.join(metadata_dir, 'fasta_rename_accession_to_gisaid_id.csv' )
    rename_fasta.to_csv(outfile, index = False)
    print('  *** WRITING RENAME FASTA FILE TO:')
    print('    %s\n' % outfile)
    
    outfile = os.path.join(metadata_dir, 'fasta_rename_accession_to_gisaid_id.tsv' )
    rename_fasta.to_csv(outfile, index = False, sep = '\t')
    print('  *** WRITING RENAME FASTA FILE TO:')
    print('    %s\n' % outfile)
    
    
    print('  FINISHED! Check warnings and output files\n\n')