__title__     = "Integrated Fragmentation Suite"
__author__    = "Swaroop Kunapuli"
__copyright__ = "Laboratoire de Modelisation et Simulations Moleculaires"
__maintainer__= "Swaroop Kunapuli"
__email__     = "svskunapuli@unistra.fr"
__version__   = "1.0"
__status__    = "Prototype"

import tempfile
import os
import sys
import numpy as np
import subprocess
import time
from Modules.IntegrationModule import FragmentationGraph
from Modules.IntegrationModule import Cluster_Combination_Object


n=int(input("Number of fundamental moieties: "))
print(n)
xyz_files = []
fundamental_moieties = []
charges = []
uhf = []

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
Global_Fragment_Database_dict={}

for number_of_precursors in range(0,2):
    number_of_moieties=[]
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
    
    res_parent = [i+j for i,j in zip([str(x) for x in parent_ion[0][:]],parent_ion[1][:])]
    name_parent = '_'.join(res_parent)

    daughter_ions, levels = FragmentationGraph(parent_ion,n,charges)
    Parent_Conformers_Database = []
    Fragments_Database = []
    if name_parent in Global_Fragment_Database_dict:
        Fragments_Database.extend(Global_Fragment_Database_dict[name_parent])
        print("Cluster already sampled!")
    if name_parent not in Global_Fragment_Database_dict:
        parent_ion_object=Cluster_Combination_Object(parent_ion,xyz_files,n)
        output_xyz, output_dir = parent_ion_object.parent_packmol_xyz_file_generation()
        parent_ion_object.crest_sampling(output_xyz=output_xyz,output_dir=output_dir,CONFORMERS_DATABASE=Parent_Conformers_Database)
        Fragments_Database.extend(Parent_Conformers_Database)
        Global_Fragment_Database_dict[name_parent]=Parent_Conformers_Database
        os.chdir("../")

    M3C_INPUT_FILE=open("{}.m3c".format(name_parent),'a+') 
    with open("M3C_TEMPLATE.m3c",'r') as M3C_TEMPLATE:
        M3C_INPUT_FILE.writelines(M3C_TEMPLATE.readlines())
    
    for ion in daughter_ions[:-1]:
        Conformers_Database=[]
        res = [i+j for i,j in zip([str(x) for x in ion[1][:]],ion[2][:])]
        name = '_'.join(res)
        if name in Global_Fragment_Database_dict:
            Fragments_Database.extend(Global_Fragment_Database_dict[name])
            print("Cluster already sampled!")        
        if name not in Global_Fragment_Database_dict: 
            ion_object=Cluster_Combination_Object(ion,xyz_files,n)
            output_xyz, output_dir = ion_object.packmol_xyz_file_generation()
            ion_object.crest_sampling(output_xyz=output_xyz,output_dir=output_dir,CONFORMERS_DATABASE=Conformers_Database)
            Fragments_Database.extend(Conformers_Database)
            Global_Fragment_Database_dict[name]=Conformers_Database
            os.chdir("../")    
    M3C_INPUT_FILE.writelines(line + '\n' for line in Fragments_Database)
    M3C_INPUT_FILE.write("END FRAGMENTS_DATABASE")
    M3C_INPUT_FILE.close()
    
