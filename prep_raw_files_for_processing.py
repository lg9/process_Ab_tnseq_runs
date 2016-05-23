# prep_raw_files_for_processing.py
# (C) Larry Gallagher 2016
#
# Written for Python 2.6

import common_parameters
import os
import re

# ------------------------------------------------------------------------------

def move_fastq_files():
    ''' Looks for folders of appropriate name (e.g, 160229-160240_S1-99) in raw_files_dir.
        If found:
         1. Copies folder to storage_dir (doesn't replace)
         2. Renames and moves them to working_dir
             e.g., 160119_S1-26/CB4_S1_L001_R1_001.fastq.gz
             becomes:  [working_dir]/160119_CB4.fastq.gz
    '''

    # Regex patterns
    dirpatt = re.compile(r'^(\d{6}(-\d+)?)_S\d')
    fqfilepatt = re.compile(r'^([a-zA-Z0-9]+)_S\d+_.*((\.fastq)|(\.fq))(\.gz)?')

    # process folders
    cwd = os.getcwd()
    os.chdir(common_parameters.raw_files_dir)
    dirlist = os.listdir('.')
    for diritem in dirlist:
        m = re.match(dirpatt, diritem)
        if (os.path.isdir(diritem) and m):
            # it is a directory with name of the expected pattern; process it.
            dirdate = None
            if m:
                dirdate = m.group(1)
            if not dirdate:
                print "Could not parse directory name, " + diritem
                continue
            # make sure raw directory and contents are backed up in storage_dir
            print "Copying directory " + diritem + " to storage location..."
            args = ['cp', '-r', '-n', diritem + '/', common_parameters.storage_dir]
            call(args)
            # rename and move run files to working_dir
            os.chdir(diritem)
            itemdirlist = os.listdir('.')
            for fqitem in itemdirlist:
                m = re.search(fqfilepatt, fqitem)
                if m:
                    sample = m.group(1)
                    extension = m.group(2)
                    if m.group(5):
                        extension += m.group(5)
                    newname = sample + '_' + dirname + extension
                    destination = working_dir + newname
                    print "   executing command: mv " \
                          + fqitem + " " + destination
                    args = ['mv', fqitem, destination]
                    call(args)
            # remove the (hopefully) empty directory
            os.chdir(common_parameters.raw_files_dir)
            args = ['rmdir', rundir + r'/']
            call(args)  # will see error if dir is not empty
            # in future, perhaps add exception handling for this call(args)
    os.chdir(cwd)

def main():
    move_fastq_files()

if __name__ == "__main__":
    main()


