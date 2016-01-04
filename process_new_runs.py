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
mapping_meta_data_file = r'/home/lg/work_Ab/mapping_meta_data.txt'
out_annot_file = r'/home/lg/work_Ab/1N_rcmp.xls'
out_tab_file = r'/home/lg/work_Ab/1N_TbGn.xls'

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
pat_annotation = r'/home/lg/work_Ab/all_features_chromosome.ptt,' + \
                 r'/home/lg/work_Ab/all_features_p1.ptt,' + \
                 r'/home/lg/work_Ab/all_features_p2.ptt,' + \
                 r'/home/lg/work_Ab/all_features_p3.ptt'
pat_options = [ '--reffile=' + pm_reference, \
                '--annofiles=' + pat_annotation, \
                '--outfile_anno=' + out_annot_file, \
                '--outfile_tab=' + out_tab_file, \
                '--workingdir=' + tnseq_wd ]

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
            args = ['rm', unzipped_file]
            call(args)
    # write meta data
    if not os.path.exists(mapping_meta_data_file):
        mdo = open(mapping_meta_data_file, 'w')
        mdo.write('Sample\tTotal reads\tTn matches\tMapped reads\tPositions\tMerged\n')
        mdo.close()
    mdo = open(mapping_meta_data_file, 'a')
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

def tabulate_samples(sample_names):
    cwd = os.getcwd()
    os.chdir(working_dir)
    print 'Annotating and tabulating samples...'
    args = ['python', tnseq_program_dir + 'process_annotate_tabulate.py'] + pat_options \
           + [sum_mg_norm_dir + s + '_trim_sum_mg_norm.txt' for s in sample_names]
    call(args)
    os.chdir(cwd)
