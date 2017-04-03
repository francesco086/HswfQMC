#!/usr/bin/env python

# Script that, given a SR_energies.dat file, find the minimum energy and correspondin standard deviation.
# Just run
#   ./find_minimum_from_SR_energies.py ottimizzazione/SR_energies.dat

# For getting the minimum energy from the bash, use
#   ./find_minimum_from_SR_energies.py ottimizzazione/SR_energies.dat | grep "minE=" | sed "s/minE= //"
#
# For getting the corresponding standard deviation, use
#   ./find_minimum_from_SR_energies.py ottimizzazione/SR_energies.dat | grep "dEmin=" | sed "s/dEmin= //"

import sys
import numpy as np

#print "find_SR_minimum_energy.py is working on ", sys.argv[1]

SRenergies = np.loadtxt(sys.argv[1], skiprows=0, usecols=(1,2))

minE=SRenergies[0][0]
dEmin=SRenergies[0][1]

for i in range (0,len(SRenergies)):
   if ((SRenergies[i][0]+1.5*SRenergies[i][1]<minE-1.5*dEmin)&(SRenergies[i][1]<dEmin*10.)):
      minE=SRenergies[i][0]
      dEmin=SRenergies[i][1]

print "minE=", minE
print "dEmin=",dEmin
