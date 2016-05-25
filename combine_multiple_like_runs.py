# comgine_multiple_like_runs.py
# (C) Larry Gallagher 2016
#
# For auto processing of A.b. tn-seq run data
# Written for Python 2.6

import os, re
from subprocess import call
from operator import itemgetter
from math import log


# NOTE:  the subroutines in this script need to be edited to work as a standalone program.


# ------------------------------------------------------------------------------

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

print 'This script needs to be edited to work as a stand-alone program.'
