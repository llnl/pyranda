#!/bin/local/bin/python
import os
import sys
import glob


# Select files for code-to-code transforms
pwd = os.path.dirname(os.path.abspath(__file__))
DIRS = []
DIRS.append(pwd + "/../../parcop")  ### Add source file dirs here ###


#### Addd your Fortran suffix support here
fsuff = ".f90"

try:
    option = int(sys.argv[1])
except Exception:
    print("Error: specify 1- forward transform, or 2-backward")
    exit()


# Option-1: Make conversion (safe to call multple times)
if option == 1:
    # Option 1- Convert files
    Files = []
    for D in DIRS:
        Files += glob.glob(D + "/*%s" % fsuff)

    # Not go through each and convert it
    for ff in Files:
        try:
            os.system("python %s/funroller.py %s %s %s" % (pwd, ff, fsuff, 0))
        except Exception:
            print("Error processing file:%s" % ff)


# Option-2: Revert space (safe to call multple times)
if option == 2:
    # Option 2- Revert files
    Files = []
    for D in DIRS:
        Files += glob.glob(D + "/*.fexl")

    for ff in Files:
        cmd = "mv -f %s %s" % (ff, ff.replace(fsuff + ".fexl", fsuff))
        os.system(cmd)
