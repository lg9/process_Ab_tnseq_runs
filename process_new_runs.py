# process_nwe_runs.py
# (C) Larry Gallagher 2016
#
# For auto processing of A.b. tn-seq run data

import os, re
from subprocess import call

# PARAMETERS and CONSTANTS -----------------------------------------------------
home_dir = r'/home/lg/'
storage_dir = r'/home/manoil-data/lg/Tn-seq_Ab/'
working_dir = r'/home/lg/work_Ab/'
tnseq_wd = r'/home/lg/work_Ab/work/'
tnseq_program_dir = r'/home/lg/Tn-seq-1.1.1/python/'
sum_mg_norm_dir = r'/home/lg/work_Ab/sum_mg_norm_files/'

# process_map
pm_reference = r'/home/lg/work_Ab/AB5075UW_allreplicons.fna'
pm_options = [  '-j', \                         # very tn end by read1
                '--tn_end=AGACAG', \            # for T26
                '--normfactor=10000000', \      # set to 10M
                '-s',                           # merge slips
                '--workingdir=' + tnseq_wd ]    # working dir for tn-seq scripts
'''
pm_options = [  '-j', '--tn_end=AGACAG', '--normfactor=10000000', '-s', '--workingdir=' + tnseq_wd ]
'''

# ------------------------------------------------------------------------------

def move_fastq_files():
    cwd = os.getcwd()
    os.chdir(home_dir)
    dirlist = os.listdir('.')
    dirpatt = re.compile(r'^\d*-?\d+_S\d+')
    fqfilepatt = re.compile(r'^\d*_?(.+)_S\d+.*((\.fastq)|(\.fq))(\.gz)?')
    for rundir in dirlist:
        if (os.path.isdir(rundir) and dirpatt.match(rundir)):
            # it is a directory with name of the expected pattern; process it.
            dirname = None
            m = re.match(r'^(\d{6}-\d+)_S\d', rundir)
            if m:
                dirname = m.group(1)
            else:
                m = re.match(r'^(\d{6})_S\d', rundir)
                if m:
                    dirname = m.group(1)
            if not dirname:
                print "Could not parse directory name, " + rundir
                continue
            print "Copying directory " + rundir + " to storage location..."
            args = ['cp', '-r', '-u', rundir + '/', storage_dir]
            call(args)
            os.chdir(rundir)
            itemdirlist = os.listdir('.')
            for item in itemdirlist:
                m = re.search(fqfilepatt, item)
                if m:
                    sample = m.group(1)
                    extension = m.group(2) + m.group(5)
                    newname = dirname + '_' + sample + extension
                    destination = working_dir + newname
                    print "   executing command: mv " \
                          + item + " " + destination
                    args = ['mv', item, destination]
                    call(args)
            os.chdir(home_dir)
            args = ['rmdir', rundir + r'/']
            call(args)  # will see error if dir is not empty
            # in future, perhaps add exception handling for this call(args)
    os.chdir(cwd)

def get_meta_data(map_logfile):
    sample_md = {}
    tot_reads_pattern = re.compile(r'Total reads processed:\s+(\d+)\s?')
    tn_match_pattern = re.compile(r'matching tn end seq:\s+(\d+)\s?')
    mapped_reads_pattern = re.compile(r'mapped reads analyzed:\s+(\d+)\s?')
    positions_pattern = re.compile(r'Positions written:\s+(\d+)\s?')
    mg_positions_pattern = re.compile(r'(\d+) positions written to.+mg\.txt')
    lfh = open(map_logfile, 'r')
    for line in lfh.readlines():
        #print line
        if 'Tot_reads' not in sample_md:
            m = tot_reads_pattern.search(line)
            if m:
                sample_md['Tot_reads'] = int(m.group(1))
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

def process_samples():
    cwd = os.getcwd()
    os.chdir(working_dir)
    dirlist = os.listdir('.')
    #sample_pattern = re.compile(r'^(\d{6}-?\d*_.+)((\.fastq)|(\.fq))(\.gz)?')
    sample_pattern = re.compile(r'^(((\d{6}-?\d*)_(.+))((\.fastq)|(\.fq)))(\.gz)?')
    meta_data = {}      # key: sample name,  value: dict of meta data
    for item in dirlist:
        # check that it's a sample file:
        m = sample_pattern.search(item)
        if m:
            unzipped_file = m.group(1)
            sample_name = m.group(2)
            date_tag = m.group(3)
            sample_tag = m.group(4)
            gzip = m.group(8)
            map_logfile = sample_name + '.map.log'
            print 'Processing sample: ', item
            # unzip if necessary
            if gzip:
                print '  unzipping...'
                args = ['gzip', '-d', item]
                call(args)
            # PROCESS the SAMPLE
            # run process_map
            #       in future: use pipe or python 2.7 check_output
            print '  running process_map.py...'
            args = ['python', tnseq_program_dir + 'process_map.py'] \
                   + ['-r', pm_reference] \
                   + pm_options + [unzipped_file]
            lfh = open(map_logfile, 'w')
            call(args, stdout=lfh)
            lfh.close()
            # capture data from the log file of process_map.py
            print '  capturing meta data...'
            meta_data[sample_name] = get_meta_data(map_logfile)
            # clean up
            print '  cleaning up files...'
            call(['rm', map_logfile])
            if not os.path.isdir(sum_mg_norm_dir):
                os.mkdir(sum_mg_norm_dir)
            os.chdir(tnseq_wd)
            args = ['cp', sample_name \
                    + '_trim_sum_mg_norm.txt', sum_mg_norm_dir]
            call(args)
            del_pattern = re.compile(r'^'+sample_name)
            files = os.listdir('.')
            for file in files:
                if del_pattern.match(file):
                    args = ['rm', file]
                    call(args)
            os.chdir(working_dir)
    os.chdir(cwd)
    # remove fastq files
    return meta_data

