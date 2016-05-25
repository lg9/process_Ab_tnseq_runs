# average_sum_files.py
# (C) Larry Gallagher 2016
#
# For auto processing of A.b. tn-seq run data
# Written for Python 2.6

import sys
import re
from operator import itemgetter
import optparse

# ------------------------------------------------------------------------------

def init_options():
    usage = "usage: %prog [options] sum_file_1 sum_file_2 ..."
    parser = optparse.OptionParser(prog=sys.argv[0],
                          usage=usage,
                          add_help_option=True)
    parser.add_option("-o", "--outsumfile", action="store", type="string", dest="outsumfile",
                      help="path to output _sum file (averages read counts from the input _sum files)")

    opts, args = parser.parse_args()
    if (opts.outsumfile is None or len(args) == 0):
        parser.print_help()
        exit(1)

    return opts, args

def average_sum_files(newsumfile, sum_files):
    '''
    Reads in the locations and read counts from multiple _sum files
    (e.g., EC3_100235_trim_sum_mg_norm.txt and EC3_100250_trim_sum_mg_norm.txt),
    and creates a new _sum file with reads per location averaged from the input files.
    '''

    # IN FUTURE,
    #  add ability to accept only locations with a minimun number of identifications
    #  add ability to work with _all.txt and _q0.txt input files

    # Regex patterns
    sumfile_patt = re.compile(r'([a-zA-Z0-9\-\+]+)_([0-9\-\+]+)' +
                              r'((_trim)?_sum(_mg)?(_norm)?((_all)|(_q0))?.txt)')
    sfpatt_grpno_type = 3

    # check that _sum files are of the same type (merged, normalized, etc.)
    sumfile_type = sumfile_patt.match(sum_files[0]).group(sfpatt_grpno_type)
    for sf in sum_files:
        sft = sumfile_patt.match(sf).group(sfpatt_grpno_type)
        if sft != sumfile_type:
            print 'All _sum files must be of same processing type (e.g., _trim_sum_mg_norm.txt)'
            return None
    # check that the new file is appropriately named
    nsft = sumfile_patt.match(newsumfile).group(sfpatt_grpno_type)
    if nsft != sumfile_type:
        print 'Specified name for new file does not match processing type (' + \
              sumfile_type + ')'
        return None

    # load all the data
    filecount = len(sum_files)
    sumfile_data = dict()   # key: (replicon, pos, dir), value:  list of (anyQ, q0) read count value tuples
    header = ''
    for sf in sum_files:
        inf = open(sf, 'r')
        header = inf.readline()
        for line in inf.readlines():
            fields = line.strip().split('\t')
            rpd = (fields[0], int(fields[1]), fields[2])
            rds = (float(fields[3]), float(fields[4]))
            if rpd not in sumfile_data:
                sumfile_data[rpd] = list
            sumfile_data[rpd].append(rds)
        inf.close()

    # make a sorted list of all the (replicons, loc, dir) tuples
    rpd_list = sumfile_data.keys()
    rpd_list.sort(key = itemgetter(2))      # first sort by the tertiary key, dir
    rpd_list.sort(key = itemgetter(1))      # next sort by secondary key, effpos
    rpd_list.sort(key = itemgetter(0))      # finally sort by primary key, replicon

    # calculate and write averages
    outf = open(newsumfile, 'w')
    outf.write(header)
    for rpd in rpd_list:
        replicon = rpd[0]
        position = str(rpd[1])
        direction = rpd[2]
        aQrdsum = 0
        q0rdsum = 0
        for rds in sumfile_data[rpd]:
            aQrdsum += rds[0]
            q0rdsum += rds[1]
        aQ_avg_rds = aQrdsum / filecount
        q0_avg_rds = q0rdsum / filecount
        # convert integer float values to integer values
        if aQ_avg_rds == int(aQ_avg_rds):
            aQ_avg_rds = int(aQ_avg_rds)
        if q0_avg_rds == int(q0_avg_rds):
            q0_avg_rds = int(q0_avg_rds)
        aQ_reads = str(aQ_avg_rds)
        q0_reads = str(q0_avg_rds)

        # write the averaged data
        lineout = '\t'.join([replicon, position, direction, aQ_reads, q0_reads])
        outf.write(lineout)

    outf.close()
    return True


def main():
    opts, args = init_options()
    average_sum_files(opts.outsumfile, args)
    
    
if __name__ == "__main__":
    main()
