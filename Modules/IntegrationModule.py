import subprocess
import os

def FragmentationGraph(parent_ion,n,charges):
    levels=0
    [number_of_moieties,fundamental_moieties,total_charge,levels]=parent_ion
    daughter_ions=[]
    seen = set()
    #For the first level of daugther ions produced from parent ion
    for i in range(0,n):
        d_number_of_moieties=number_of_moieties[:]
        d_number_of_moieties[i]=d_number_of_moieties[i]-1
        daughter_ions.append([number_of_moieties,d_number_of_moieties,fundamental_moieties,sum([x*y for x,y in zip(d_number_of_moieties,charges)]),levels+1])
        del(d_number_of_moieties)
    levels=levels+1
    #For the second level of daughter ions produced from those daughter ions (unique listing)
    max_levels=50
    for i in range(2,max_levels):
        for item in daughter_ions:
            if item[4]==i-1:
                for j in range(0,n):
                    d_number_of_moieties=item[1][:]
                    d_number_of_moieties[j]=d_number_of_moieties[j]-1
                    if d_number_of_moieties[j]>=0:
                        k=tuple(d_number_of_moieties)
                        if k not in seen:
                            seen.add(k)
                            daughter_ions.append([item[1][:],d_number_of_moieties,fundamental_moieties,sum([x*y for x,y in zip(d_number_of_moieties,charges)]),i])
                    del(d_number_of_moieties)
    number_of_daughter_ions=0
    for item in daughter_ions:
        print(item)
        number_of_daughter_ions=number_of_daughter_ions+1
    print(number_of_daughter_ions)

# using CREST software to sample the conformational space
def crest_sampling(parent_ion,solute_xyz,solvent_xyz,number_of_solvent_molecules,xyz_files):
    [number_of_moieties,fundamental_moieties,levels]=parent_ion
    number_of_moieties_str=[]
    for item in number_of_moieties:
        number_of_moieties_str.append(str(item))
    if os.path.exists(fundamental_moieties[1]):
        subprocess.run(["cp",xyz_files[1],fundamental_moieties[1]])
        subprocess.run(["cp",xyz_files[0],fundamental_moieties[1]])
    else:
        subprocess.run(["mkdir",fundamental_moieties[1]])
        subprocess.run(["cp",xyz_files[0],fundamental_moieties[1]])
        subprocess.run(["cp",xyz_files[1],fundamental_moieties[1]])
    os.chdir(fundamental_moieties[1])

    subprocess.run(["crest",solute_xyz,"--qcg",solvent_xyz,"--nsolv",number_of_solvent_molecules,"--xtbiff","/home/software/xtbiff","--gfnff","--T","5","--nofix"])
    #os.system("cd grow")
    print("NCI + HESS (1) OR NCI_THEN_HESS (2)?")
    nci_hess_decision=input()
    if nci_hess_decision=="1":
        os.system("cp grow/cluster_optimized.xyz .")
        subprocess.run(["crest","cluster_optimized.xyz","--chrg", str(int(charges[1])+(int(charges[0])*int(number_of_moieties[0]))) ,"--nci","--prop","hess","--T","5", "&&"],capture_output=True)
        GetFrequencies('PROP')
    if nci_hess_decision=="2":
        os.system("cp grow/cluster_optimized.xyz .")
        subprocess.run(["crest","cluster_optimized.xyz","--chrg",str(int(charges[1])+(int(charges[0])*int(number_of_moieties[0]))),"--nci","--T","5","&&"])
        subprocess.run(["crest","crest_best.xyz","--chrg", str(int(charges[1])+(int(charges[0])*int(number_of_moieties[0]))),"--for","crest_conformers.xyz","--prop","hess","--T","5","&&"],capture_output=True)
        GetFrequencies('PROP')

# Calculate the frequencies and the total energy of the system and output the 1) RXYZ files of each conformer 
def GetFrequencies(properties_directory):
    f = open("cre_members",'r',)
    lines=f.readlines()
    number_of_conformers=int(lines[0])
    ##subprocess.run(["cd",properties_directory])
    #Printinf rxyz files for conformers from frequencies, electronic energies and coordinates in PROP/
    for i in range(1,number_of_conformers+1):
        RXYZ_FILE=open(f'conformer{i}.rxyz','a+')
        struc_xyz=open(f'{properties_directory}/TMPCONF{i}/struc.xyz','r')
        RXYZ_FILE.write(struc_xyz.read())
        with open(f'{properties_directory}/TMPCONF{i}/vibspectrum','r') as f:
            frequencies=[]
            lines=f.readlines()
            l=0
            for line in lines:
                words=line.split()
                if len(words)>=2 and words[1]=='a':
                    frequencies.append(words[2])
                    l=l+1
            RXYZ_FILE.write('\n')
            RXYZ_FILE.write("FREQUENCIES %d \n" %l) 
            RXYZ_FILE.write('\n'.join(str(freq) for freq in frequencies))
        RXYZ_FILE.close()
                
# Append to M3C Fragment Database in the format prescribed
def MakeM3CFragmentDatabase():
    M3C_Fragmentation_Database=open('M3C_Fragmentation_Database.txt','a+')
    
