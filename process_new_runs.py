# process_nwe_runs.py
# (C) Larry Gallagher 2016
#
# For auto processing of A.b. tn-seq run data

import os, re
from subprocess import call

# PARAMETERS and CONSTANTS -----------------------------------------------------
storage_dir = r'/home/manoil-data/lg/Tn-seq_Ab/'
working_dir = r'~/work_Ab/'
# ------------------------------------------------------------------------------

def move_fastq_files():
    dirlist = os.listdir('.')
    dirpatt = re.compile(r'^\d*-?\d+_S\d+')
    fqfilepatt = re.compile(r'^\d*_?(.+)_S\d+.*((\.fastq)|(\.fq))(\.gz)?')
    for item in dirlist:
        if (os.path.isdir(item) and dirpatt.match(item)):
            # it is a directory of the expected name; process it.
            print "Copying directory " + item + " to storage location..."
            args = ['cp', '-r', item + '/', storage_dir]
            # exception handle the next command to cover if already copied
            #call(args)
            dirname = None
            m = re.match(r'^(\d{6}-\d+)_S\d', item)
            if m:
                dirname = m.group(1)
            else:
                m = re.match(r'^(\d{6})_S\d', item)
                if m:
                    dirname = m.group(1)
            if not dirname:
                print("Could not parse directory", item)
                continue
            os.chdir(item)
            itemdirlist = os.listdir('.')
            for item in itemdirlist:
                m = re.search(fqfilepatt, item)
                if m:
                    sample = m.group(1)
                    extension = m.group(2) + m.group(5)
                    newname = dirname + '_' + sample + extension
                
                    print "mv " + item + " " + working_dir + newname
                    
            os.chdir("..")
            
