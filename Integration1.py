import os
import sys
import numpy as np
import subprocess
import Modules

#os.system("crest --help")
#subprocess.run(["crest","--help"])

n=int(input("Number of fundamental moieties: "))
print(n)
xyz_files = []
fundamental_moieties = []
charges = []
uhf = []
number_of_moieties = []
for i in range(1,n+1):
    print("If .xyz files of the Fundamental Moieties are in this directory, enter their name. If not, enter the location + name of the files")
    xyz_files.append(input())
    print("Give it a label:")
    fundamental_moieties.append(input())
    print("Its charge:")
    charges.append(int(input()))
    print("Its Multiplicity:")
    uhf.append(input())
print(fundamental_moieties)
print("Enter the composition of the precursor ion")
for i in range(1,n+1):
    print(fundamental_moieties[i-1],":")
    number_of_moieties.append(int(input()))
print("Parent non-covalent Ion at level 0 of fragmentation has following molecular entities : ")
for i in range(0,n):
   print(number_of_moieties[i],fundamental_moieties[i])
levels=0
parent_ion=[number_of_moieties,fundamental_moieties, sum([x*y for x,y in zip(number_of_moieties,charges)]),levels]
print(parent_ion)

Modules.FragmentationGraph(parent_ion,n,charges)
#print("")
#Modules.crest_sampling(parent_ion,xyz_files[1],xyz_files[0],"1",xyz_files)
