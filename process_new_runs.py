# process_nwe_runs.py
# (C) Larry Gallagher 2016
#
# For auto processing of A.b. tn-seq run data
# Written for Python 2.6

import os, re
from subprocess import call
from operator import itemgetter
from math import log

# PARAMETERS and CONSTANTS -----------------------------------------------------
home_dir = r'/home/lg/'
storage_dir = r'/home/manoil-data/lg/Tn-seq_Ab/'
working_dir = r'/home/lg/work_Ab2/'
tnseq_wd = r'/home/lg/work_Ab2/work/'
tnseq_program_dir = r'/home/lg/Tn-seq-1.1.2/python/'
sum_mg_norm_dir = r'/home/lg/work_Ab2/sum_mg_norm_files/'
mapping_meta_data_file = r'/home/lg/work_Ab2/mapping_meta_data.txt'
out_annot_file = r'/home/lg/work_Ab2/1N_rcmp.xls'
out_tab_file = r'/home/lg/work_Ab2/1N_TbGn.xls'
wig_file_dir = r'/home/lg/work_Ab2/wig_files/'

# process_map
pm_reference = r'/home/lg/work_Ab/AB5075UW_allreplicons.fna'
pm_options = [  '-j', \                         # very tn end by read1
                '--tn_end=AGACAG', \            # for T26
                '--normfactor=10000000', \      # set to 10M
                '-s',                           # merge slips
                '-k',                           # sequencing is from 'back end'
                '--workingdir=' + tnseq_wd ]    # working dir for tn-seq scripts
'''
pm_options = [  '-j', '--tn_end=AGACAG', '--normfactor=10000000', '-s', '-k', '--workingdir=' + tnseq_wd ]
'''

# process_annotate_tabulate
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
    ''' Looks for directories like '160101_S1-2' or '151228-160101_S1-4', etc.
        in home directory, then:
        1.  copies them to the storage_dir
        2.  renames and moves them to the working_dir
            e.g., 160119_S1-26/CB4_S1_L001_R1_001.fastq.gz
            becomes:  [working_dir]/160119_CB4.fastq.gz
    '''
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
    ''' For all run files in working_dir (e.g., 160119_CC4.fastq.gz):
         1. unzips if necessary
         2. runs process_map.py, capturing and saving output meta data
         3. copies _trim_sum_mg_norm.txt file to sum_mg_norm_dir
         4. removes all working files, including source fastq file
        Then, adds meta data for each run to mapping_meta_data_file
    '''
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

def sample_info_f_runfile_name(filename):
    ''' Sub function for extracting sample info from a run filename.
        Called by combine_multiple_runs()
    '''
    fn_pattern = re.compile(r'(^\d*_?(.+)_S(\d+).*((\.fastq)|(\.fq)))(\.gz)?')
    m = fn_pattern.match(filename)
    if not m:
        return None
    sample = m.group(2)
    info = { sample: {} }
    info[sample]['unzipped_file'] = m.group(1)
    info[sample]['S_no'] = int(m.group(3))
    info[sample]['fq_type'] = m.group(4)
    info[sample]['zipped'] = m.group(7)
    return info

def combine_multiple_runs(folder_list, new_folder, directory=None):
    '''
    Function to concatenate run fastq files
    from multiple runs (in separate folders) of the same samples.
    Combined runs are created in the specified new_folder.
    
    Does this in the home directory or the specified one
    '''
    cwd = os.getcwd()
    if not directory:
        directory = home_dir
    os.chdir(directory)
    for folder in folder_list:
        if not os.path.isdir(folder):
            print "source directory, " + folder + " not found."
            return None
    if not os.path.isdir(new_folder):
        os.mkdir(new_folder)
    # create lists of files to concatenate per sample
    print "Reading and unzipping contents of source folders..."
    sample_cat_lists = {}
    sample_S_nos = {}
    for folder in folder_list:
        filelist = os.listdir(folder)
        for fn in filelist:
            rni = sample_info_f_runfile_name(fn)
            # skip if it's not a sample run file
            if not rni:     
                continue
            sample = rni.keys()[0]
            if sample not in sample_cat_lists:
                sample_cat_lists[sample] = []
                sample_S_nos[sample] = []
            sample_path = folder + '/' + rni[sample]['unzipped_file']
            sample_cat_lists[sample].append(sample_path)
            if rni[sample]['S_no'] not in sample_S_nos[sample]:
                sample_S_nos[sample].append(rni[sample]['S_no'])
            if rni[sample]['zipped']:
                print '   unzipping ' + folder + '/' + fn
                args = ['gzip', '-d', folder + '/' + fn]
                call(args)
    # concatenate the files for each sample and place in new_folder
    print "Concatenating and (re-)zipping files..."
    for sample in sample_cat_lists.keys():
        print '   concatenating files for sample ' + sample
        s_no = 'S' + str(sample_S_nos[sample][0])
        if len(sample_S_nos[sample]) > 1:
            for i in range(1,len(sample_S_nos[sample])):
                s_no += '-' + str(sample_S_nos[sample][i])
        newfilename = sample + '_' + s_no + '.fastq'
        destfileandpath = new_folder + '/' + newfilename
        fo = open(destfileandpath, 'w')
        args = ['cat'] + sample_cat_lists[sample]
        call(args, stdout=fo)
        fo.close()
        files_to_zip = sample_cat_lists[sample] + [destfileandpath]
        for f in files_to_zip:
            print '      zipping ' + f
            args = ['gzip', f]
            call(args)
    print "Done"
    os.chdir(cwd)

def tabulate_samples(sample_names):
    ''' Runs process_annotate_tabulate.py
        on the corresponding _trim_sum_mg_norm.txt files
        from a list of specified samples.
    '''
    cwd = os.getcwd()
    os.chdir(working_dir)
    print 'Annotating and tabulating samples...'
    args = ['python', tnseq_program_dir + 'process_annotate_tabulate.py'] \
           + pat_options + [sum_mg_norm_dir + s + '_trim_sum_mg_norm.txt' \
                            for s in sample_names]
    call(args)
    os.chdir(cwd)

def make_wig_files(rcmpfile=out_annot_file, makeqn0=False, logtfm=False):
    '''
    Makes wig files of the runs in a rcmp file.
    wig files placed in wig_file_dir.
    Options:
        makeqn0: set to T to make wig files for qn0 reads
        logtfm: set to T to log2 transform data
    '''
    # Read data from rcmp file and place in data structure
    print "   reading data..."
    inf = open(rcmpfile, 'r')
    header = inf.readline()
    h_flds = header.split('\t')
    h_fld_count = len(h_flds)
    rplcn_fld_no3 = None
    effpos_fld_no = None
    dir_fld_no = None
    notes_fld_no = None
    for i in range(h_fld_count):
        if 'eplicon' in h_flds[i]:
            rplcn_fld_no = i
        if 'ffective' in h_flds[i]:
            effpos_fld_no = i
        if 'Dir' in h_flds[i]:
            dir_fld_no = i
        if 'Notes' in h_flds[i]:
            notes_fld_no = i
    samples = list()
    sample_count = 0
    wigfile_names = list()
    sample_fld_nos = list()
    sample_pattern = re.compile(\
        r'(\d+-?\d*_[a-zA-Z0-9]+)(_trim)?_sum(_mg)?(_norm)?((_all)|(_q0))?.txt')
    for i in range(notes_fld_no + 1, h_fld_count):
        m = sample_pattern.match(h_flds[i])
        sample = m.group(1)
        # check to see if it's a q0 sample
        if sample in samples and m.group(7):
            # if so, record the number of samples and quit the loop
            sample_count = len(samples)
            break
        samples.append(sample)
        # determine the wig file name
        wigfile_name = sample + '_'
        if m.group(3):
            wigfile_name += 'mg'
        if m.group(4):
            wigfile_name += '1N'
        if makeqn0:
            wigfile_name += 'qn0'
        else:
            wigfile_name += 'aQ'
        if logtfm:
            wigfile_name += '_log2'
        wigfile_name += '.wig'
        wigfile_names.append(wigfile_name)
        sample_fld_nos.append(i)
    # read data into data structure:
    #   a list, each entry a list of replicon, pos, dir and sample read counts
    data = list()
    min_non0_val = float(10000000.0)
    for line in inf.readlines():
        flds = line.split('\t')
        datum = list()
        datum.append(flds[rplcn_fld_no])
        datum.append(int(flds[effpos_fld_no]))
        datum.append(flds[dir_fld_no])
        for i in sample_fld_nos:
            rds = float(flds[i])
            if rds > 0 and rds < min_non0_val:
                min_non0_val = rds
            if makeqn0:
                rds -= float(flds[i + sample_count])
            datum.append(rds)
        data.append(datum)
    inf.close()
    print "      lines read: " + str(len(data))
    # if logtfm set to T, first re-normalize all read counts so that the minimum
    # non-zero value is 2, then log2-transform all non-zero values.
    if logtfm:
        print "   log-transforming data..."
        r = range(3, len(data[0]))
        for d in range(len(data)):
            for e in r:
                if data[d][e] > 0:
                    data[d][e] = log( (data[d][e] * 2 / min_non0_val), 2 )
    # sort the data
    data.sort(key = itemgetter(2))      # first sort by the tertiary key, dir
    data.sort(key = itemgetter(1))      # next sort by secondary key, effpos
    data.sort(key = itemgetter(0))      # finally sort by primary key, replicon
    # Make .wig files
    print "   writing wig files..."
    # create wig_file_dir if it doesn't exist
    if not os.path.isdir(wig_file_dir):
        os.mkdir(wig_file_dir)
    # create array of filehandles for the various wig files
    fhs = list()
    for wfn in wigfile_names:
        wfpn = wig_file_dir + wfn
        fh = open(wfpn, 'w')
        fhs.append(fh)
    # write the data
    r = range(sample_count)
    chmsm = ''
    for datum in data:
        if datum[0] != chmsm:
            chmsm = datum[0]
            new_replicon_line = 'variableStep\tchrom=' + chmsm + '\n'
            for fh in fhs:
                fh.write(new_replicon_line)
        effpos = datum[1]
        if datum[2].strip() == 'F':
            dir_factor = 1
        else:
            dir_factor = -1
        for sn in r:    # sample number
            if datum[sn+3] > 0:
                fhs[sn].write(str(effpos) + '\t' \
                              + str(dir_factor * datum[sn+3]) + '\n')
    # close filehandles
    for fh in fhs:
        fh.close()
