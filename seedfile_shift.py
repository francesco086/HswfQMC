#!/usr/bin/env python

from subprocess import call

seedfin = open('randomseed1.d', 'r')
seedfout = open('randomseed_new.d', 'w')

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

call(["mv","randomseed_new.d","randomseed1.d"])
