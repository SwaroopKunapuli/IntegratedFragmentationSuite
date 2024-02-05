import Modules.pypackmol as pyp
import os
import subprocess

def FragmentationGraph(parent_ion,n,charges):
    levels=0
    [number_of_moieties,fundamental_moieties,total_charge,levels]=parent_ion
    global daughter_ions
    daughter_ions=[]
    seen = set()
    #For the first level of daugther ions produced from parent ion
    for i in range(0,n):
        d_number_of_moieties=number_of_moieties[:]
        d_number_of_moieties[i]=d_number_of_moieties[i]-1
        daughter_ions.append([[number_of_moieties],d_number_of_moieties,fundamental_moieties,sum([x*y for x,y in zip(d_number_of_moieties,charges)]),levels+1])
        del(d_number_of_moieties)
    levels=levels+1
    break_out_flag=False
    while break_out_flag==False: 
        for item in daughter_ions:
            if item[4]==levels: 
                for j in range(0,n):
                    d_number_of_moieties=item[1][:]
                    d_number_of_moieties[j]=d_number_of_moieties[j]-1
                    if all(numbers_d>=0 for numbers_d in d_number_of_moieties):
                        k=tuple(d_number_of_moieties)
                        if k not in seen:
                            seen.add(k)
                            daughter_ions.append([[item[1][:]],d_number_of_moieties,fundamental_moieties,sum([x*y for x,y in zip(d_number_of_moieties,charges)]),levels+1])
                        else:
                            for index, value in enumerate(daughter_ions):
                                if value[1][:]==d_number_of_moieties:
                                        daughter_ions[index][0].append(item[1][:])
                    if all(numbers_d<=0 for numbers_d in d_number_of_moieties):
                        break_out_flag=True
                    del(d_number_of_moieties)
        levels=levels+1
    print("Daughter ions (cluster combinations) from this parent ion are: ", *daughter_ions, sep='\n')
    print("Total number of ions (cluster combinations): " + str(len(daughter_ions)))
    print("Number of levels in the Fragmentation Graph : " + str(levels))
    return daughter_ions, levels


                
# Class to generate the following information related to each parent/daughter ion from the fragmentation graph
# INPUT:: 1) Object is the the parent/daugther ion 
#         2) Packmol generation of the .xyz file of the parent/daughter ion with pypackmol interdace
# PROCESSING: 4) IMPLEMENTING THE FUNCTION FOR a) (Optional) Optimisation of the structure
#                                              b) NCI conformer sampling 
#                                              c) Hessian Calculation for FREQUENCIES AND TOTAL ENERGY 
# RETURNING: 5) RETURNING the a) .rxyz files 
#                             b) append to M3C Fragment Database 
#                             c) best cluster conformer to go to their parent ion object
    
class Cluster_Combination_Object(object):
    ion_dict={}
    def __init__(self,ion,xyz_files,number_of_fundamental_moieties):
        """
        Generating packed .xyz file from Packmol with pypackmol wrapper for each of the ions in the parent/daughter ions list.
        
        """
        self.ion = ion
        name_ion="{}{}_{}{}_{}{}".format(self.ion[1][0],self.ion[2][0],self.ion[1][1],self.ion[2][1],self.ion[1][2],self.ion[2][2])
        Cluster_Combination_Object.ion_dict[name_ion] = self
        self.xyz_files=xyz_files
        self.number_of_fundamental_moieties =number_of_fundamental_moieties

    def packmol_xyz_file_generation(self):
        pm = pyp.Packmol(dimension=10)
        for i in range(0,self.number_of_fundamental_moieties):
            if self.ion[1][i]>0:
                pm.add_structure(self.xyz_files[i],count=self.ion[1][i],input_format="xyz")
        global output_xyz
        global output_dir
        print(self.ion[1][0],self.ion[2][0],self.ion[1][1],self.ion[2][1],self.ion[1][2],self.ion[2][2])
        output_dir="{}{}_{}{}_{}{}".format(self.ion[1][0],self.ion[2][0],self.ion[1][1],self.ion[2][1],self.ion[1][2],self.ion[2][2])
        output_xyz="{}{}_{}{}_{}{}.xyz".format(self.ion[1][0],self.ion[2][0],self.ion[1][1],self.ion[2][1],self.ion[1][2],self.ion[2][2])
        result=pm.pack(output=str(output_xyz))
        return output_xyz, output_dir
    
    # Calculate the frequencies and the total energy of the system and output the 1) RXYZ files of each conformer 
    def GetFrequencies(self,properties_directory):
        if os.path.exists('cre_members'):
            f = open("cre_members",'r')
            lines=f.readlines()
            number_of_conformers=int(lines[0])
            Negative_conformer_list_file=open("Negative_Frequency_conformers_{}.txt".format(output_dir),'a+')
            #subprocess.run(["cd",properties_directory])
            #Printinf rxyz files for conformers from frequencies, electronic energies and coordinates in PROP/
            for i in range(1,number_of_conformers+1):
                RXYZ_FILE=open(f'{output_dir}_{i}.rxyz','a+')
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
                            if float(words[2])<0:
                                Negative_conformer_list_file.write("Negative frequencies in conformer{}".format(i))
                    RXYZ_FILE.write('\n')
                    RXYZ_FILE.write("FREQUENCIES %d \n" %l) 
                    RXYZ_FILE.write('\n'.join(str(freq) for freq in frequencies))
                RXYZ_FILE.close()
            Negative_conformer_list_file.close()
        else:
            print("No cre_members file!!")

    # using CREST software to sample the conformational space
    def crest_sampling(self,output_xyz,output_dir ):
        if os.path.exists(output_dir):
            subprocess.run(["cp",output_xyz,output_dir])
        else:
            subprocess.run(["mkdir",output_dir])
            subprocess.run(["cp",output_xyz,output_dir])
        os.chdir(output_dir)
        crest_nci_output=open("crest_nci_output.txt","w")
        subprocess.run(["crest",output_xyz,"--chrg",str(self.ion[3]),"--nci","--gfn2//gfnff","--T","5","--mdlen","0.5","&&"],stdout=crest_nci_output)
        crest_hess_output=open('crest_hess_output.txt','w')
        subprocess.run(["crest",output_xyz,"--chrg",str(self.ion[3]),"--gfn2//gfnff","--for","crest_conformers.xyz","--prop","hess","--T","5","&&"],stdout=crest_hess_output)
        if os.path.exists('PROP'):
            self.GetFrequencies('PROP')
    #def add_to_fragment_database(self,):