#! /usr/bin/env python

import sys 
import subprocess
import shlex
import pandas as pd

def download_fastq(run, line_dat):
    # get extension	
    # this is sloppy. can't think of a better way right now. 
    print(line_dat.iloc[0])
    # check if file exists
    ps = subprocess.Popen(shlex.split(f'gsutil -m ls {line_dat["fastq_dir"]}'), 
        stdout=subprocess.PIPE)
    ps_head = subprocess.check_output(shlex.split('head -n 1'), stdin=ps.stdout)
    ext = ps_head.decode('utf-8').strip().split('.')[-1]
    if ext == 'gz':
        # concatenated file name
        outfile = f'{run}_fastq/CO-CDPHE-{line_dat.iloc[0]}.fastq.gz'
        # concatenate barcoded reads to single file
        with open(outfile, 'w') as fp:
            subprocess.run(shlex.split(f'gsutil -m cat {line_dat["fastq_dir"]}/*.fastq.gz'), stdout=fp)
    elif ext=='fastq':
        # concatenated file name
        outfile = f'{run}_fastq/CO-CDPHE-{line_dat.iloc[0]}.fastq.gz'
        # concatenate barcoded reads to single file
        with open(outfile, 'w') as fp:
            subprocess.run(shlex.split(f'gsutil -m cat {line_dat["fastq_dir"]}/*.fastq'), stdout=fp)
#     else:
#         raise Exception('either .fastq.gz or .fastq file expected')


if __name__ == "__main__":
    for run in sys.argv[2:]:
        print(run)
        # read in filtered results file
        df = pd.read_csv(sys.argv[1], sep = '\t', dtype = {'accession_id' : object})
        crit = df.seq_run == sys.argv[2]
        df = df[crit]
        df = df.reset_index(drop = True)
        
        include = \
            set(df.iloc[:,0].astype(str))

        # makes fastq directory if it doesn't exist
        subprocess.run(shlex.split(f'mkdir -p {run}_fastq'))
        
        # download terra data table
        subprocess.run(shlex.split(f'gsutil -m cp gs://covid_terra/{run}/{run}_terra_data_table.tsv .'))
        
        # read in terra data table
        dat = pd.read_csv(f'{run}_terra_data_table.tsv', sep='\t', dtype = {0: object})
        
        # subset 
        dat = dat[dat.iloc[:,0].isin(include)]
        print('len of filtered terra datatable = %d rows' % dat.shape[0])
        dat.apply(lambda k: download_fastq(run, k), axis=1)


