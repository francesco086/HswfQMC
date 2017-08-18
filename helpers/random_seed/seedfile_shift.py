#!/usr/bin/env python
# This script takes a random seed file for HswfQMC as argument
# and moves the first seed line to the end of the file 
#
# Usage: seedfile_shift.py <randomseed.d>

from subprocess import call
import sys

oldfile = sys.argv[1]
nfsplit = oldfile.split('.')

if nfsplit[-1] == 'd' and len(nfsplit)>1:
    nfsplit[-1] = 'new'
else:
    raise ValueError('Random seed filename is expected to end with .d')

newfile = ''.join(nfsplit)

seedfin = open(oldfile, 'r')
seedfout = open(newfile, 'w')

lcount = 0
for line in seedfin:
    lcount += 1
    if lcount == 2:
        swapline = line
    else:
        seedfout.write(line)

seedfout.write(swapline)
seedfout.close()
seedfin.close()

call(['mv', newfile, oldfile])
