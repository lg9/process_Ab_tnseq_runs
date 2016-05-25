# process_nwe_runs.py
# (C) Larry Gallagher 2016
#
# For auto processing of A.b. tn-seq run data
# Written for Python 2.6

import common_parameters
import os, re
from subprocess import call

# ------------------------------------------------------------------------------

def get_meta_data(map_logfile):
    ''' Sub function to get meta data from a saved map_logfile.
        Called by map_samples()
    '''
    sample_md = {}
    tot_reads_pattern = re.compile(r'Total reads processed:\s+(\d+)\s?')
    tn_match_pattern = re.compile(r'matching tn end seq:\s+(\d+)\s?')
    mapped_reads_pattern = re.compile(r'mapped reads analyzed:\s+(\d+)\s?')
    positions_pattern = re.compile(r'Positions written:\s+(\d+)\s?')
    mg_positions_pattern = re.compile(r'(\d+) positions written to.+mg\.txt')
    lfh = open(map_logfile, 'r')
    for line in lfh.readlines():
        #print line
        if 'Total_reads' not in sample_md:
            m = tot_reads_pattern.search(line)
            if m:
                sample_md['Total_reads'] = int(m.group(1))
        if 'Tn_matches' not in sample_md:
            m = tn_match_pattern.search(line)
            if m:
                sample_md['Tn_matches'] = int(m.group(1))
        if 'Mapped_reads' not in sample_md:
            m = mapped_reads_pattern.search(line)
            if m:
                sample_md['Mapped_reads'] = int(m.group(1))
        if 'Positions' not in sample_md:
            m = positions_pattern.search(line)
            if m:
                sample_md['Positions'] = int(m.group(1))
        if 'Merged_positions' not in sample_md:
            m = mg_positions_pattern.search(line)
            if m:
                sample_md['Merged_positions'] = int(m.group(1))
    lfh.close()
    return sample_md

def map_samples():
    ''' For all run files in working_dir (e.g., CC4_160119.fastq.gz):
         1. unzip if necessary
         2. runs process_map.py, capturing and saving output meta data
         3. copies _trim_sum_mg_norm.txt file to sum_mg_norm_dir
         4. removes all working files, including source fastq file
        Then, adds meta data for each run to mapping_meta_data_file
    '''
    
    # Regex patterns
    sample_pattern = re.compile(r'((([a-zA-Z0-9\-\+]+)_(\d{6}[0-9\-\+]*))((\.fastq)|(\.fq)))(\.gz)?')

    # process samples
    cwd = os.getcwd()
    os.chdir(common_parameters.working_dir)
    dirlist = os.listdir('.')
    meta_data = {}      # key: sample name,  value: dict of meta data
    for item in dirlist:
        # check that it's a sample file:
        m = sample_pattern.search(item)
        if m:
            # parse info from name
            unzipped_file = m.group(1)
            sample_name = m.group(2)
            sample_tag = m.group(3)
            date_tag = m.group(4)
            gzipped = m.group(8)
            map_logfile = sample_name + '.map.log'
            print 'Processing sample: ', item
            # unzip if necessary
            if gzipped:
                print '  unzipping...'
                args = ['gzip', '-d', item]
                call(args)
            # PROCESS the SAMPLE
            # run process_map
            #       in future: use pipe or python 2.7 check_output
            print '  running process_map.py...'
            args = ['python', common_parameters.tnseq_program_dir + 'process_map.py'] \
                   + ['-r', common_parameters.reference_genome] \
                   + common_parameters.pm_options + [unzipped_file]
            log_fh = open(map_logfile, 'w')
            call(args, stdout=log_fh)
            log_fh.close()
            # capture data from the log file of process_map.py
            print '  capturing meta data...'
            meta_data[sample_name] = get_meta_data(map_logfile)
            # clean up
            print '  cleaning up files...'
            call(['rm', map_logfile])
            if not os.path.isdir(common_parameters.sum_mg_norm_dir):
                os.mkdir(common_parameters.sum_mg_norm_dir)
            os.chdir(common_parameters.tnseq_wd)
            args = ['cp', sample_name \
                    + '_trim_sum_mg_norm.txt', common_parameters.sum_mg_norm_dir]
            call(args)
            del_pattern = re.compile(r'^'+sample_name)
            files = os.listdir('.')
            for file in files:
                if del_pattern.match(file):
                    args = ['rm', file]
                    call(args)
            os.chdir(common_parameters.working_dir)
            args = ['rm', unzipped_file]
            call(args)
    # write meta data
    if not os.path.exists(common_parameters.mapping_meta_data_file):
        mdo = open(common_parameters.mapping_meta_data_file, 'w')
        mdo.write('Sample\tTotal reads\tTn matches\tMapped reads\tPositions\tMerged\n')
        mdo.close()
    mdo = open(common_parameters.mapping_meta_data_file, 'a')
    for sample in sorted(meta_data.keys()):
        mdo.write('\t'.join([ sample, \
                             str(meta_data[sample]['Total_reads']), \
                             str(meta_data[sample]['Tn_matches']), \
                             str(meta_data[sample]['Mapped_reads']), \
                             str(meta_data[sample]['Positions']), \
                             str(meta_data[sample]['Merged_positions']) ] ) + '\n' )
    mdo.close()
    os.chdir(cwd)
    return meta_data

def main():
    md = map_samples()

if __name__ == "__main__":
    main()

