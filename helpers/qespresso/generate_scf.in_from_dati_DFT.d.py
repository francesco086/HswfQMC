#! /usr/bin/env python

import subprocess
import sys
import time
import re

print "Start the scf.in generator ..."

### GET THE INFORMATION NECESSARY TO BUILD THE INPUT FILE FOR QUANTUM ESPRESSO

#Get the path to the folder where the pseudo potential to be used is stored
mystring=subprocess.check_output("   cat ../../dati_DFT.d | grep 'PSEUDO_DIR=' | sed 's/PSEUDO_DIR=//'   ",shell=True)
mystring=re.sub("[\"'\n]","",mystring)
#print "The path to the pseudo potential folder is " mystring
pseudo_dir=mystring

#Get pseudo potential to be used
mystring=subprocess.check_output("   cat ../../dati_DFT.d | grep 'PSEUDO=' | sed 's/PSEUDO=//'   ",shell=True)
mystring=re.sub("[\"'\n]","",mystring)
#print "The pseudo potential to use is " mystring
pseudo=mystring

#Get the cut off energy
mystring=subprocess.check_output("   cat ../../dati_DFT.d | grep 'ECUTWFC=' | sed 's/ECUTWFC=//'   ",shell=True)
mystring=re.sub("[\"'\n]","",mystring)
ecutwfc=mystring

#Get the K-points to use
mystring=subprocess.check_output("   cat ../../dati_DFT.d | grep 'K_POINTS=' | sed 's/K_POINTS=//'   ",shell=True)
mystring=re.sub("[\"'\n]","",mystring)
k_points=mystring

#Get the K-grid to use
mystring=subprocess.check_output("   cat ../../dati_DFT.d | grep 'K_GRID=' | sed 's/K_GRID=//'   ",shell=True)
mystring=re.sub("[\"'\n]","",mystring)
k_grid=mystring

#Get the size of the box from the temproray file L.d
mystring=subprocess.check_output("cat L.d",shell=True)
mylist=mystring.split()
L=[eval(x) for x in mylist]
#print "L = ", L

#Get the crystal structure from the temporary file crystal.d
mystring=subprocess.check_output("cat crystal.d",shell=True)
mylist=mystring.split()
foo=[eval(x) for x in mylist]
i=0
r_crystal=[]
while i<len(foo):
    r_crystal.append(foo[i:i+3])
    i+=3
#print "r_crystal = ", r_crystal

#Set the number of atoms
N=len(r_crystal)


### WRITE THE INPUT FILE FOR QUANTUM ESPRESSO

#Open the file scf.in, which will be the input file for pw.x (quantum espresso)
QE_input_file=open('scf.in', 'w')

#Write on QE_input_file the CONTROL namelist
mystring="&CONTROL\n"
mystring+="calculation = \"scf\"\n"
mystring+="restart_mode = \"from_scratch\"\n"
mystring+="pseudo_dir = \""+pseudo_dir+"\"\n"
mystring+="prefix = \"OUT\"\n"
mystring+="wf_collect = .TRUE.\n"
mystring+="/\n"
QE_input_file.write(mystring)

#Write on QE_input_file the SYSTEM namelist
mystring="&SYSTEM\n"
mystring+="ibrav = 0\n"
#mystring+="celldm(1) = 1.0\n"
mystring+="nat  = "+str(N)+"\n"
mystring+="ntyp  =  1\n"
mystring+="ecutwfc  =   "+ecutwfc+"\n"
mystring+="occupations  =  \"smearing\"\n"
mystring+="smearing = \"gaussian\"\n"
mystring+="degauss = 0.005\n"
mystring+="/\n"
QE_input_file.write(mystring)

#Write on QE_input_file the ELECTRONS namelist
mystring="&ELECTRONS\n"
mystring+="conv_thr = 1.0d-6\n"
mystring+="mixing_beta = 0.5\n"
mystring+="/\n"
QE_input_file.write(mystring)

#Write on QE_input_file the CELL_PARAMETERS part
mystring="CELL_PARAMETERS (cubic) bohr\n"
mystring+=str(L[0])+"   0.0000000   0.0000000\n"
mystring+="0.0000000   "+str(L[1])+"   0.0000000\n"
mystring+="0.0000000   0.0000000   "+str(L[2])+"\n"
QE_input_file.write(mystring)

#Write on QE_input_file the ATOMIC_SPECIES part
mystring="ATOMIC_SPECIES\n"
mystring+="H  1.00794 H.coulomb-ae.UPF \n"
QE_input_file.write(mystring)

#Write on QE_input_file the ATOMIC_POSITIONS
mystring="ATOMIC_POSITIONS (crystal)\n"
#Append the atomic positions
i=0
while i<N:
    line_to_write="H"
    j=0
    while j<len(r_crystal[i]):
        line_to_write+="   "
        line_to_write+=str(r_crystal[i][j])
        j+=1
    #print line_to_write
    line_to_write+="\n"
    mystring+=line_to_write
    i+=1
QE_input_file.write(mystring)

#Write on QE_input_file the K_POINTS
mystring="K_POINTS "+k_points+"\n"
mystring+=k_grid+" 0 0 0\n"
QE_input_file.write(mystring)

print "... scf.in generated succesfully!"


