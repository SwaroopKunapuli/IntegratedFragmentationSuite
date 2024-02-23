import Modules.pypackmol as pyp
import os
import subprocess
import time

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
    __slots__ = ['ion','xyz_files','number_of_fundamental_moieties'] 
    def __init__(self,ion,xyz_files,number_of_fundamental_moieties):
        """
        Generating packed .xyz file from Packmol with pypackmol wrapper for each of the ions in the parent/daughter ions list.
        
        """
        self.ion = ion
        self.xyz_files=xyz_files
        self.number_of_fundamental_moieties =number_of_fundamental_moieties
    
    def parent_packmol_xyz_file_generation(self):
        pm = pyp.Packmol(dimension=10)
        for i in range(0,self.number_of_fundamental_moieties):
            if self.ion[0][i]>0:
                pm.add_structure(self.xyz_files[i],count=self.ion[0][i],input_format="xyz")
        global output_xyz
        global output_dir
        res=[i+j for i,j in zip([str(x) for x in self.ion[0][:]],self.ion[1][:])]
        print(res)
        output_dir='_'.join(res)
        output_xyz="{}.xyz".format(output_dir)
        result=pm.pack(output=str(output_xyz))
        return output_xyz, output_dir
    
    def packmol_xyz_file_generation(self):
        pm = pyp.Packmol(dimension=10)
        for i in range(0,self.number_of_fundamental_moieties):
            if self.ion[1][i]>0:
                pm.add_structure(self.xyz_files[i],count=self.ion[1][i],input_format="xyz")
        global output_xyz
        global output_dir
        res=[i+j for i,j in zip([str(x) for x in self.ion[1][:]],self.ion[2][:])]
        print(res)
        output_dir='_'.join(res)
        output_xyz="{}.xyz".format(output_dir)
        result=pm.pack(output=str(output_xyz))
        return output_xyz, output_dir
    
    def AutoMeKin_llcals(self,list_xyz,i,name,charge):
        subprocess.run(["mkdir","{}_{}_AutoMeKin".format(name,i)])
        subprocess.run(["cp","AutoMekin_Template.dat","{}_{}_AutoMeKin".format(name,i)])
        subprocess.run(["cp","qcore_template","{}_{}_AutoMeKin".format(name,i)])
        os.chdir("{}_{}_AutoMeKin".format(name,i))
        AutoMeKin_Input_File=open("{}_{}.dat".format(name,i),'a+')
        AutoMeKin_Output_File=open("output.txt",'a+')
        with open("{}_{}.xyz".format(name,i),'a+') as XYZ_FILE:
            XYZ_FILE.writelines(list_xyz)
        with open("AutoMekin_Template.dat",'r') as AutoMeKin_Template:
            AutoMeKin_Template_lines=AutoMeKin_Template.readlines()
            AutoMeKin_Input_File.writelines(AutoMeKin_Template_lines[0])
            AutoMeKin_Input_File.write("molecule {}_{} \n".format(name,i))
            AutoMeKin_Input_File.write("charge {} \n".format(charge))
            AutoMeKin_Input_File.writelines(AutoMeKin_Template_lines[2:])
        AutoMeKin_Input_File.close()
        subprocess_llcalcs_procees = subprocess.Popen("nohup llcalcs.sh {}_{}.dat 5 5 5 >llcals.log 2>&1 &".format(name,i),stdout=subprocess.PIPE,shell=True)
        (subprocess_llcalcs_output, subprocess_llcalcs_error) = subprocess_llcalcs_procees.communicate()
        subprocess_llcalcs_status = subprocess_llcalcs_procees.wait()
        print("Subprocess output : {}".format(subprocess_llcalcs_output))
        TS_DE=None
        file_path="FINAL_LL_{}_{}/TSinfo".format(name,i)
        while not os.path.exists(file_path):
            time.sleep(200)
        if os.path.isfile(file_path):
            with open(file_path,'r') as TSinfo:
                TSinfo_readlines=TSinfo.readlines()
                for line in list(TSinfo_readlines):
                    words=line.split()
                    if len(words)>=2 and words[0]=='1':
                        TS_DE=float(words[1])*0.0433641
        else:
            raise ValueError("No transition state file ()".format(file_path))
        os.chdir("../")
        AutoMeKin_Output_File.close()
        return TS_DE       



    # Calculate the frequencies and the total energy of the system and output the 1) RXYZ files of each conformer 
    def GetFrequencies(self,properties_directory,CONFORMERS_DATABASE,charge):
        if os.path.exists('cre_members'):
            f = open("cre_members",'r')
            lines=f.readlines()
            number_of_conformers=int(lines[0])
            Negative_conformer_list_file=open("Negative_Frequency_conformers_{}.txt".format(output_dir),'a+')
            for i in range(1,number_of_conformers+1):
                Negative_conformer=False
                RXYZ_FILE=open(f'{output_dir}_{i}.rxyz','a+')
                if os.path.exists('{}/TMPCONF{}/xtbopt.xyz'.format(properties_directory,i)):
                    xtb_opt_xyz=open(f'{properties_directory}/TMPCONF{i}/xtbopt.xyz','r')
                    xtb_opt_xyz_lines=xtb_opt_xyz.readlines()
                    RXYZ_FILE.write(xtb_opt_xyz_lines[0])
                    TS_DE = self.AutoMeKin_llcals(xtb_opt_xyz_lines,i,output_dir,charge)
                    xtb_out=open(f'{properties_directory}/TMPCONF{i}/xtb.out','r')
                    xtb_out_lines=xtb_out.readlines()
                    for line in reversed(list(xtb_out_lines)): 
                        words=line.split()
                        if len(words)>3 and words[1]=='total' and words[2]=='energy':
                            energy = words[3]
                    RXYZ_FILE.write("Energy = {}\n".format(energy))
                    RXYZ_FILE.write("".join(line for line in xtb_opt_xyz_lines[2:]))
                if os.path.exists('{}/TMPCONF{}/vibspectrum'.format(properties_directory,i)):
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
                                    Negative_conformer=True
                        if Negative_conformer==True:
                            Negative_conformer_list_file.write("Negative frequencies in conformer{} \n".format(i))
                        if Negative_conformer==False:
                            CONFORMERS_DATABASE.append("{}_{}  {} 1 0 1 {}/{}_{}.rxyz {} {}".format(output_dir,i,self.ion[3],output_dir,output_dir,i,str(27.2114*float(energy)),TS_DE))
                        RXYZ_FILE.write('\n')
                        RXYZ_FILE.write("FREQUENCIES %d \n" %l) 
                        RXYZ_FILE.write('\n'.join(str(freq) for freq in frequencies))
                elif os.path.exists('{}/TMPCONF{}/NOT_CONVERGED'.format(properties_directory,i)):
                    Negative_conformer_list_file.write("NOT CONVERGED IN conformer{} \n".format(i))
                else:
                    print("Some error!")
                RXYZ_FILE.close()
            Negative_conformer_list_file.close()
        else:
            print("No cre_members file!!")
            
    # using CREST software to sample the conformational space
    def crest_sampling(self,output_xyz,output_dir,CONFORMERS_DATABASE):
        if os.path.exists(output_dir):
            subprocess.run(["cp",output_xyz,output_dir])
            subprocess.run(["cp","AutoMekin_Template.dat",output_dir])
            subprocess.run(["cp","qcore_template",output_dir])
        else:
            subprocess.run(["mkdir",output_dir])
            subprocess.run(["cp",output_xyz,output_dir])
            subprocess.run(["cp","AutoMekin_Template.dat",output_dir])
            subprocess.run(["cp","qcore_template",output_dir])
        os.chdir(output_dir)
        crest_nci_output=open("crest_nci_output.txt","w")
        subprocess.run(["crest",output_xyz,"--chrg",str(self.ion[3]),"--nci","--gfn2//gfnff","--T","5","--mdlen","0.5","&&"],stdout=crest_nci_output)
        crest_hess_output=open('crest_hess_output.txt','w')
        subprocess.run(["crest",output_xyz,"--chrg",str(self.ion[3]),"--for","crest_conformers.xyz","--prop","ohess","--gfn2","--T","5","&&"],stdout=crest_hess_output)
        if os.path.exists('PROP'):
            self.GetFrequencies('PROP',CONFORMERS_DATABASE,self.ion[3])

    def parent_crest_sampling(self,output_xyz,output_dir,CONFORMERS_DATABASE):
        if os.path.exists(output_dir):
            subprocess.run(["cp",output_xyz,output_dir])
            subprocess.run(["cp","AutoMekin_Template.dat",output_dir])
            subprocess.run(["cp","qcore_template",output_dir])
        else:
            subprocess.run(["mkdir",output_dir])
            subprocess.run(["cp",output_xyz,output_dir])
            subprocess.run(["cp","AutoMekin_Template.dat",output_dir])
            subprocess.run(["cp","qcore_template",output_dir])
        os.chdir(output_dir)
        crest_nci_output=open("crest_nci_output.txt","w")
        subprocess.run(["crest",output_xyz,"--chrg",str(self.ion[2]),"--nci","--gfn2//gfnff","--T","5","--mdlen","0.5","&&"],stdout=crest_nci_output)
        crest_hess_output=open('crest_hess_output.txt','w')
        subprocess.run(["crest",output_xyz,"--chrg",str(self.ion[2]),"--for","crest_conformers.xyz","--prop","ohess","--gfn2","--T","5","&&"],stdout=crest_hess_output)
        if os.path.exists('PROP'):
            self.GetFrequencies('PROP',CONFORMERS_DATABASE,self.ion[2])
        
    