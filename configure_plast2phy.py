#!/usr/bin/env python

# Checks requirements for plast2phy

# DPENDENCIES
# Python modules
#   sh (http://amoffat.github.io/sh/)
#   Biopython (python-biopython)
# Programs/scripts in PATH or ./bin/
#   mafft (apt-get install mafft) [http://mafft.cbrc.jp/alignment/software/linux.html]
#   trimal v.1.2 () [http://trimal.cgenomics.org/]
#   raxml v.8.x.x() [https://github.com/stamatak/standard-RAxML]

import imp
import os
import sys


def check_modules(modules):
    '''Checks if the nessecary modules are installed'''
    for module in modules:
        try:
            imp.find_module(module)
        except ImportError:
            if module == 'Bio':
                sys.exit("ERROR: This pipeline requires Biopython. In Ubuntu run " +\
                    "sudo apt-get install python-biopython. If that doesn't work " +\
                    "see Biopython's website for further install methods http://biopython.org/")
            elif module == 'sh':
                sys.exit("ERROR: This pipeline requires the sh library (http://amoffat.github.io/sh/)" +\
                    "the easiest way of installig it is through pip, pip install sh " +\
                    "In Ubuntu the pip package is under python-pip")


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def main():
    '''Checks for required programs and modules'''
    modules = ['Bio', 'sh']
    required_programs = ['mafft', 'raxmlHPC']
    check_modules(modules)
    for program in required_programs:
        path_true = which(program)
        if path_true == None:
            sys.exit("%s is missing, please make sure %s is in PATH" % (program, program))
    print "All dependencies are met, the pipeline should work."


main()