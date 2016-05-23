# common_parameters.py
# (C) Larry Gallagher 2016
#
# Written for Python 2.6


# FILE and DIRECTORY names and locations -----------------------------------------------------
home_dir = r'/home/lg/'
raw_files_dir = home_dir
storage_dir = r'/home/manoil-data/lg/Tn-seq_Ab/'
working_dir = r'/home/lg/work_Ab2/'
tnseq_program_dir = r'/home/lg/Tn-seq-1.1.2/python/'
tnseq_wd = r'/home/lg/work_Ab2/work/'
sum_mg_norm_dir = r'/home/lg/work_Ab2/sum_mg_norm_files/'
mapping_meta_data_file = r'/home/lg/work_Ab2/mapping_meta_data.txt'
out_annot_file = r'/home/lg/work_Ab2/1N_rcmp.xls'
out_tab_file = r'/home/lg/work_Ab2/1N_TbGn.xls'

wig_file_dir = r'/home/lg/work_Ab2/wig_files/'

# process_map
reference_genome = working_dir + r'AB5075UW_allreplicons.fna'
pm_options = [  '-j', \                         # very tn end by read1
                '--tn_end=AGACAG', \            # for T26
                '--normfactor=10000000', \      # set to 10M
                '-s',                           # merge slips
                '-k',                           # sequencing is from 'back end'
                '--workingdir=' + tnseq_wd ]    # working dir for tn-seq scripts

# process_annotate_tabulate
pat_annotation = working_dir + r'all_features_chromosome.ptt,' + \
                 working_dir + r'all_features_p1.ptt,' + \
                 working_dir + r'all_features_p2.ptt,' + \
                 working_dir + r'all_features_p3.ptt'
pat_options = [ '--reffile=' + reference_genome, \
                '--annofiles=' + pat_annotation, \
                '--outfile_anno=' + out_annot_file, \
                '--outfile_tab=' + out_tab_file, \
                '--workingdir=' + tnseq_wd ]

# ------------------------------------------------------------------------------

