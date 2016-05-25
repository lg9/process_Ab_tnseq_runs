# make_wig_files.py
# (C) Larry Gallagher 2016
#
# For auto processing of A.b. tn-seq run data
# Written for Python 2.6

import common_parameters
import os
import re
from operator import itemgetter
from math import log

# ------------------------------------------------------------------------------

def make_wig_files(rcmpfile=common_parameters.out_annot_file, makeqn0=False, logtfm=False):
    '''
    Makes wig files of the runs in a rcmp file.
    wig files placed in wig_file_dir.
    Options:
        makeqn0: set to T to make wig files for qn0 reads
        logtfm: set to T to log2 transform data
    '''

    # Regex
    sample_pattern = re.compile(\
        r'([a-zA-Z0-9\-\+]+_\d+([\-\+]\d+)?)(_trim)?_sum(_mg)?(_norm)?((_all)|(_q0))?.txt')

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
    for i in range(notes_fld_no + 1, h_fld_count):
        m = sample_pattern.match(h_flds[i])
        sample = m.group(1)
        # check to see if it's a q0 sample
        if sample in samples and m.group(8):
            # if so, record the number of samples and quit the loop
            sample_count = len(samples)
            break
        samples.append(sample)
        # determine the wig file name
        wigfile_name = sample + '_'
        if m.group(4):
            wigfile_name += 'mg'
        if m.group(5):
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
    if not os.path.isdir(common_parameters.wig_file_dir):
        os.mkdir(common_parameters.wig_file_dir)
    # create array of filehandles for the various wig files
    fhs = list()
    for wfn in wigfile_names:
        wfpn = common_parameters.wig_file_dir + wfn
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

def main():
    print 'making linear .wig files...'
    make_wig_files()
    print 'making log2-transformed .wig files...'
    make_wig_files(logtfm=True)
    print 'done.'
    
if __name__ == "__main__":
    main()
