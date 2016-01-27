'''
make_wig_files.py
(C) 2016 Larry Gallagher

Written for python 2.7

For making wig files from 1N_rcmp.xls file
'''



#import sys
import argparse
import re
from math import log
from operator import itemgetter


# PARAMETERS and CONSTANTS -----------------------------------------------------

DEF_RCMP_FILE = '1N_rcmp.xls'

#------------------------------------------------------------------------------

def init_options():
    parser = argparse.ArgumentParser( \
        description='Create wig files from rcomp txt file - ' + \
        'for viewing hits and reads with IGV.')
    parser.add_argument("-i", "--infile", dest="infile", required=True, \
                        help='path to input file (default='+DEF_RCMP_FILE+').' \
                        default=DEF_RCMP_FILE)
#    parser.add_argument("-o", "--outfilebase", dest="outfilebase", \
#                        help="path and base name of output file (e.g. '141017_2N_').")
    parser.add_argument("-l", "--logtfm", dest="logtfm", action='store_true', \
                        help="normalize to a global non-zero min of 2 and log2-transform.")
    # in future:  add argument for capturing q0 or qn0 data
    args = parser.parse_args()
    return args

def is_even(x):
    return (x % 2) == 0

def read_rcmp_file(args):
    '''
    Read all data from rcmp file and place in data structure
    '''
    inf = open(args.infile)
    # get sample names and field numbers from header line
    header = inf.readline()
    h_fields = header.split('\t')
    h_field_count = len(h_fields)
    rpcln_f_no = None
    effpos_f_no = None
    dir_f_no = None
    notes_f_no = None
    for i in range(h_field_count):
        if 'eplicon' in h_fields[i]:
            rpcln_f_no = i
        if 'ffective' in h_fields[i]:
            effpos_f_no = i
        if 'Dir' in h_fields[i]:
            dir_f_no = i
        if 'Notes' in h_fields[i]:
            notes_f_no = i
    sample_fns = list()     # list of sample wig file names
    sample_fld_nos = list()   # correspoding list of inf field numbers
    sample_pattern = re.compile(\
        r'(\d+-?\d*_[a-zA-Z0-9]+)(_trim)?_sum(_mg)?(_norm)?((_all)|(_q0))?.txt')
    for i in range(notes_f_no + 1, h_field_count):
        m = sample_pattern.match(h_fields[i])
        if not m:       # skip if it doesn't match the sample name pattern
            continue
        sample_fn = m.group(1) + '_'
        if m.group(3):
            sample_fn += 'mg'
        if m.group(4):
            sample_fn += '1N'
        if m.group(6):
            sample_fn += 'aQ'
        if m.group(7):
            sample_fn += 'q0'
        sample_fn += '.wig'
        sample_fns.append(sample_fn)
        sample_fld_nos.append(i)
    # create header line for data
    data_wh = list()        # data with header
    data_header = ['Replicon', 'EffPos', 'Dir']
    data_header += sample_fns
    data_wh.append(data_header)
    # get read count data
    data = list()
    for line in inf.readlines():
        flds = line.split('\t')
        datum = list()
        datum.append(flds[rpcln_f_no])
        datum.append(int(flds[effpos_f_no]))
        datum.append(flds[dir_f_no])
        for i in sample_fld_nos:
            datum.append(float(flds[i]))
        data.append(datum)
    inf.close()
    # sort data

    data_wh += data
    return data_wh

# KEEP WORKING FROM HERE

def log2_like_transform_data():
    '''
    To transform data for display as log2-transformed, first re-normalize all data so that the minimum non-zero
    value is 2, then log-2 transform all the non-zero values while values of zero transform to zero.  Thus,
    the minimum read count will appear with a transformed value of 1 (log2 of 2), twice that number of reads
    will have a transformed value of 2 (log2 of 4), etc.
    '''
    global data
    # first determine overall non-zero minimum
    minval = 10000000.0
    for datum in data:
        for val in datum[3:]:
            if val > 0:
                if val < minval:
                    minval = val
    # now re-normalize and log2-transform non-zero values
    r = range(3, len(data[0]))
    for d in range(len(data)):
        for e in r:
            if data[d][e] > 0:
                data[d][e] = log( (data[d][e] * 2 / minval), 2 )

def make_wig_files(args):
    global data, sample_names
    # sort the data by chromosome, effpos and dir
    data.sort(key = itemgetter(2))      # sort by tertiary key, dir
    data.sort(key = itemgetter(1))      # now sort by secondary key, effpos
    data.sort(key = itemgetter(0))      # finally sort by primary key, replicon
    # create array of filehandles for accessing each data file; open the files for writing
    ofb = args.outfilebase
    if args.logtfm:
        ofb += '_log2_'
    fhs = list()
    for sample in sample_names:
        ofn = ofb.strip() + sample.strip() + '.wig'
        fh = open(ofn, 'w')
        fhs.append(fh)
    # parse and write the data
    r = range(len(fhs))
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
        for sn in r:
            if datum[sn+3] > 0:
                fhs[sn].write(str(effpos) + '\t' + str(dir_factor * datum[sn+3]) + '\n')
    # close filehandles
    for fh in fhs:
        fh.close()

def main():
    # Get command line options
    args = init_options()
    # Read data
    data = read_rcmp_file(args)
    if args.logtfm:
        log2_like_transform_data()
    make_wig_files(args)
 
    print("Done.")

if __name__ == "__main__":
    main()
