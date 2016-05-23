# annotate_tabulate_mapped_samples.py
# (C) Larry Gallagher 2016
#
# For auto processing of A.b. tn-seq run data
# Written for Python 2.6

import common_parameters
import os
import optparse
from subprocess import call

# ------------------------------------------------------------------------------

def init_options():
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(prog=sys.argv[0],
                                   usage=usage,
                                   add_help_option=True)
    parser.add_option("-s", "--samplenamesfile", action="store",
                      type="string", dest="samplenamesfile", required=True,
                      help="path to text file listing samples to compare")
    opts, args = parser.parse_args()

    return opts, args

def get_sample_names_f_samplenamesfile(sample_names_file):
    samples = list()
    with snf as open(sample_names_file, 'r'):
        for sl in snf.readlines():
            samples.append(sl.trim())
    return samples

def tabulate_samples(sample_names):
    ''' Runs process_annotate_tabulate.py
        on the corresponding _trim_sum_mg_norm.txt files
        from a list of specified samples.
    '''
    cwd = os.getcwd()
    os.chdir(common_parameters.working_dir)
    print 'Annotating and tabulating samples...'
    args = ['python', common_parameters.tnseq_program_dir + \
            'process_annotate_tabulate.py'] \
           + common_parameters.pat_options + \
           [common_parameters.sum_mg_norm_dir + s + \
            '_trim_sum_mg_norm.txt' for s in sample_names]
    call(args)
    os.chdir(cwd)

def main():
    opts, args = init_options()
    sample_names = get_sample_names_f_samplenamesfile(opts.samplenamesfile)
    tabulate_samples(sample_names)

if __name__ == "__main__":
    main()

