from ase import Atoms
from ase import io
from ase.build import bulk, sort
from ase.io import read
from xespresso import Espresso
from xespresso.hpxespresso import HpXEspresso
from itertools import chain
import numpy as np
import copy
import os, shutil
import sys

def update_UVscsd(U_dict,Hubbard_site):
    """
    Updates the labels of the Hubbard sites in the structure on the base of the U value.
    If there are more than 9 different Hubbard sites, letters and not numbers will be used to
    label the site becausepw.x cannot distinguish for example "Ti1" from "Ti12".
    The different number of Hubbard sites cannot exceed 26 (letters of the English alphabet)
    
    WARNING: By default pw.x can deal only up to 10 different Hubbard atoms. If you need more you need to
    change the source code and re-compile your pw.x code.
    
    U_dict: dict
        dictionary obtained with the read_Hubbard_parameters function
    Hubbard_site: str
        A str indicatint the type of the Hubbard specie
        
    Returns:
        - U: dict
            A copy of U_dict in which the new_label is modified according on the U values 
    """
    
    U = copy.deepcopy(U_dict)
    import string as STR
    alphabet = list(STR.ascii_lowercase)


    new_labels = []
    new_labels = []
    Hubbard_U = {}


    count_sites_spin0 = []
    count_sites_spin1 = []
    count_sites_spinm1 = []
    for site, info in U.items():
        if info["spin"] == "1":
            count = info["new_label"][len(Hubbard_site):]
            if count == "":
                count = 1
            else:
                count = int(count)
            count_sites_spin1.append(count)
        elif info["spin"] == "-1":
            count = info["new_label"][len(Hubbard_site):]
            if count == "":
                count = 1
            else:
                count = int(count)
            count_sites_spinm1.append(count)
        elif info["spin"] == "0":
            count = info["new_label"][len(Hubbard_site):]
            if count == "":
                count = 1
            else:
                count = int(count)
            
            count_sites_spin0.append(count)
        else:
            pass
            

    count_sites_spin0 = len(count_sites_spin0) 
    count_sites_spin1 = len(count_sites_spin1)
    count_sites_spinm1 = len(count_sites_spinm1)


    new_labels = []
    for site, info in U.items():
        if info["spin"] == "1":
            count = info["new_label"][len(Hubbard_site):]
            if count == "":
                count = 1
            else:
                count = int(count)
            new_labels.append(count)
            info["new_new_labels"] = count
 
        elif info["spin"] == "-1":
            count = info["new_label"][len(Hubbard_site):]
            if count == "":
                count = 1
            else:
                count = int(count)
            new_count = count + count_sites_spin1
            new_labels.append(new_count)
            info["new_new_labels"] = new_count

        else:
            count = info["new_label"][len(Hubbard_site):]
            if count == "":
                count = 1
            else:
                count = int(count)
            new_count = count + count_sites_spin1  + count_sites_spinm1
            new_labels.append(new_count)
            info["new_new_labels"] = new_count

            
    new_labels.sort()
    new_labels =set(new_labels)
    indeces = {}
    
    if len(new_labels) >  len(alphabet):
        sys.exit("I can differentiate maximum 26 sites")
        
    for n,i in enumerate(new_labels):
        indeces[i] = n+1
    #print(indeces)
    
    
    for site, info in U.items():
        count = info["new_new_labels"]
        indx = indeces[count]
        #print(count,indx)
        if indx < 10:
            info["new_label"] = Hubbard_site + str(indx)
        else:
            info["new_label"] = Hubbard_site + alphabet[indx]
            
        
    return U 

def update_Hubbard_parameters(U_dict,pseudopotentials,Hubbard_site,atoms): 
    """
    Updating the input_ntyp, pseudopotentials, and structure on the base of the new_labels/U values
    
    U_dict: dict
        dictionary obtained with the read_Hubbard_parameters function
        
    pseudopotentials: dict
        A dict specifing the pseudopotential name per atom type
    
    Hubbard_site: str
        A str indicatint the type of the Hubbard specie
        
    atoms: ase.atoms
        An ase structure object
        
    Returns:
        - dict with the updated structure, pseudopotentials, and input_ntyp dictionary  
    """
    Hubbard_U = dict()
    starting_magnetization = dict()
    for site, info in U_dict.items():
        label = info["new_label"]
        if label not in Hubbard_U:
            Hubbard_U[label] = float(info["U"])
            starting_magnetization[label] = float(info["spin"])
            
    remove_species = []
    for elem, pseudo in pseudopotentials.items():
        if Hubbard_site in elem:
            pseudo_name = pseudopotentials[elem]
            remove_species.append(elem)
            
    for elem in remove_species:    
        del pseudopotentials[elem]
            
    for site in Hubbard_U.keys():
        pseudopotentials[site] = pseudo_name
        
    atoms.info['species'] = atoms.get_chemical_symbols()
    
    for site, info in U_dict.items():
        atoms.info['species'][int(site)-1] = info["new_label"]
    
    return {"Hubbard_U": Hubbard_U,
            "starting_magnetization":starting_magnetization,
            "pseudopotentials"  : pseudopotentials,
           "atoms": atoms}

class Uscsd():
    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='xespresso', atoms=None, package = 'pw', parallel =  '',
                 queue = None,queue_hp = False, parallel_hp = '',max_iter= 1, iteration=1, sc_type =None,
                 insulator = None, magnetic = None, nbnd =  None, relax_type = "vc-relax",
                 parallelization_hp = False, step = None, Hubbard_site = None, conv_Uscsd = 0.01,
                 ethr_relax = None, ethr_scf = None, ethr_scf2 = None,
                 **kwargs):
        """
        This class allows to automate the self-consistent (site-dependent) calculation of U/V values
        Automating the necessary relax-scf-hp procedure. In particular, one iteration requires:
        
            - For a metal:
                1) relaxation (optional)
                2) scf calculation with "smaering"(on the optimized structure if step 1) is performed)
                3) hp calculation
                
            - For a non-magnetic insulator:
                1) relaxation (optional)
                2) scf calculation with "fixed"(on the optimized structure if step 1) is performed)
                3) hp calculation
                
            - For a magnetic insulator:
                1) relaxation (optional)
                2) scf calculation with "smaering"(on the optimized structure if step 1) is performed)
                3) scf calculation with "fixed" reading the total_magnetization and the number of bands 
                   from step 1
                4) hp calculation
                
        When the U/V values are obtained, the subsequent iteration is run using the new U/V values and eventually 
        the structure optimized in step 1) until a maximum numer of iterations (max_iter) is performed or until
        U/V convergence is reached.
        
        The class is able to detect the procedure to apply. WARNING: When sc_type = "scf" and iteration = 1,
        you need to provide values for "insulator", "magnetic", "nbnd".
        
        The class is able to restart an interrupted self-consistent calculation of U/V if for one iteration
        only the "relax" or "scf" step perfomed. 
        WARNING: restart from an incomplete "relax" step is not yet possible.
        
        
        Check:
        - "Hubbard parameters from density-functional perturbation theory", 
        Phys. Rev. B 98, 085127 (2018); arXiv:1805.01805 using ase and XEspresso, 
        of which pXEspresso is a subclass.
        - C. Ricca, I. Timrov, M. Cococcioni, N. Marzari, and U. Aschauer,
          "Self-consistent site-dependent DFT+U study of stoichiometric and defective SrMnO3",
          Phys. Rev. B 99, 094102 (2019); arXiv:1811.10858
        - C. Ricca, I. Timrov, M. Cococcioni, N. Marzari, and U. Aschauer,
          "Self-consistent DFT+U+V study of oxygen vacancies in SrTiO3",
          Phys. Rev. Research 2, 023313 (2020); arXiv:2004.04142
        
        input_data, pseudopotentials, kspacing, kpts, koffset
            Please have a look at Espresso module in ASE
            
        max_iter: int
            An integer indicating the maximum number of iteration
            
        Hubbard_site: str
            A str indicating the atomic species of the Hubbard atom for which the self-consistent (site-dependent)
            procedure is to be performd. Example: Hubbard_site = "Ti".
            
        conv_Uscsd: float
            A float indicating he threshold to be used to verify convergence of the U/V values
        
        !!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!
        Do not forget to provide in the input_data dictionry and two elements:
            "lda_plus_U" : True
            "input_ntyp": { "starting_magnetization" : {Hubbard_site: XX},
                            "Hubbard_U" : {Hubbard_site: XX}}
        """
        
        self.label_ref = label

        self.directory = os.path.split(label)[0]
        self.prefix = os.path.split(label)[1]
        self.directory_ref = self.directory
        
        self.queue_hp = queue_hp
        self.parallel_hp = parallel_hp
        self.queue = queue
        self.parallel = parallel
            
        self.step = None
        self.iteration = iteration
        self.max_iter = max_iter
        
        if sc_type == None:
            self.sc_type = "relax"
        else:
            self.sc_type = sc_type
        

        self.relax_type = relax_type

        if self.sc_type != None and self.sc_type not in ["relax", "scf"]:
            sys.exit("Please provide a str for the sc_type variable: relax or scf")

        if self.relax_type != None and self.relax_type not in ["relax", "vc-relax"]:
            sys.exit("Please provide a str for the relax_type variable: relax or vc-relax")
        if self.sc_type != None and self.sc_type not in ["relax", "scf"]:
            sys.exit("Please provide a bool for the sc_type variable: relax or scf")

        self.insulator = insulator
        self.magnetic = magnetic
        self.nbnd = nbnd

        if self.insulator != None and bool(self.insulator) == False:
            sys.exit("Please provide a bool for the insulator variable")

        if self.magnetic != None and self.magnetic not in [1,2]:
            sys.exit("Please provide an int for the magnetic  variable: 1 for spin unpolirzed and 2 for spin poolariz")

        if self.nbnd != None and type(self.nbnd) != int:
            sys.exit("Please provide an int for the nbnd variable")


        self.parallelization_hp = parallelization_hp
    
        
        self.converged = False
        
        self.kwargs = kwargs
        
        self.atoms = atoms
        
        self.UVscsd = {}
        self.etot = {}
        self.opt_structures = {}
        
        if Hubbard_site == None:
            sys.exit("You need to specify the Hubbard_site ")
        else:
            self.Hubbard_site = Hubbard_site
            
        self.conv_Uscsd = conv_Uscsd
        
        self.ethr_relax = ethr_relax
        self.ethr_scf = ethr_scf
        self.ethr_scf2 = ethr_scf2
        
        
    def converged_U(self):
        old_U = []
        new_U = []

        for site, info in self.UVscsd[self.iteration-1].items():
            old_U.append(float(info["U"]))

        for site, info in self.UVscsd[self.iteration].items():
            new_U.append(float(info["U"]))

        #comparison = np.isclose(np.array(old_U),np.array(new_U),conv_Uscsd)
        self.converged = np.allclose(np.array(old_U),np.array(new_U),self.conv_Uscsd)       
        
    def check_iteration(self):
        """
        In case previous relax-scf-hp iterations have been previously performed, this function
        identifies the last iteration and which is the last step at which that iteration was stopped
        
        WARNING:
        - Not yet able to identify if structural relaxation did not convrge and needs restart
        
        TO DO:
            - make it smarter being able to recognize only "relax_*", "scf_*", "hp_*"
        
        """
        
        folders = [name for name in os.listdir(self.directory) if os.path.isdir(os.path.join(self.directory,name))]
        for folder_path in folders:
            if len(os.listdir(os.path.join(self.directory,folder_path))) == 0: # Check is empty..
                shutil.rmtree(os.path.join(self.directory,folder_path)) # Delete..
                folders.remove(folder_path)
  
        Folders = []
        for folder_path in folders:
            if "_"  in folder_path:
                Folders.append(folder_path)
            else:
                pass
                #folders.remove(folder_path
        
        if not folders:
            dirs = []
        else:
            dirs =  [int(name.split('_')[1]) for name in Folders]
 
        if not dirs:
            self.iteration = 1
            self.step = "re-start"
        else:
            self.iteration = max(dirs)

            if os.path.exists(os.path.join(self.directory,"hp_"+str(self.iteration))):
                self.step = "hp"
            elif os.path.exists(os.path.join(self.directory,"scf_"+str(self.iteration))):
                self.step = "scf"
                print("""WARNING: For magnetic insulators, restart of the second scf step with 'fixed'
                from the previous scf step with 'smearing' is not yet implemented. In this case, the
                hp step will fail if you do not remove also the scf_$iteration folder.""")
            elif os.path.exists(os.path.join(self.directory,"relax_"+str(self.iteration))):
                self.step = "relax"
            else:
                self.step = "re-start"
                

            
    def run_sc(self):
        """
        Main function wunning the self-consistence procedure
        """
        
        if os.path.exists(self.directory) == False:
            os.mkdir(self.directory)
            print("Creating working directory")
            
        
        self.check_iteration()
        
        print("-------------------------------------------")            
        print("Initial Iteration:",self.iteration,self.step)
        print("-------------------------------------------")
        
        calc = HpXEspresso(label=self.label_ref,
                           parallel = self.parallel,
                           parallel_hp = self.parallel_hp,
                           queue = self.queue,
                           queue_hp = self.queue_hp,
                           max_iter = self.max_iter,
                           iteration = self.iteration,
                           sc_type = self.sc_type,
                           relax_type = self.relax_type,
                           insulator = self.insulator,
                           magnetic = self.magnetic,
                           nbnd = self.nbnd,
                           parallelization_hp = self.parallelization_hp,
                           step = self.step,
                           ethr_relax = self.ethr_relax,
                           ethr_scf = self.ethr_scf,
                           ethr_scf2 = self.ethr_scf2,
                           **self.kwargs)
        
        if self.iteration > 1:
            U = calc.read_Hubbard_parameters(self.iteration-1)
            self.UVscsd[self.iteration-1] =  U
            
            if self.step == 'hp':
                U = calc.read_Hubbard_parameters(self.iteration)
                self.UVscsd[self.iteration] =  U
                self.converged_U()
                print("U/V values converged: {}".format(self.converged))
                if self.converged == True:
                    if self.sc_type == "relax":
                        directory = os.path.join(self.directory_ref,"relax_"+str(self.iteration))
                    else:
                        directory = os.path.join(self.directory_ref,"scf_"+str(self.iteration))
                    calc.read_results_fromdirectory(directory)
                    self.etot[self.iteration] = calc.get_potential_energy()
                    self.opt_structures[self.iteration] = calc.atoms
                    sys.exit("Stopping: Convergence achieved!")
            
            new_U = update_UVscsd(U,self.Hubbard_site)
            new_data = update_Hubbard_parameters(new_U,self.kwargs['pseudopotentials'],self.Hubbard_site,self.atoms)
            self.atoms = new_data["atoms"]
            self.kwargs['pseudopotentials'] = new_data['pseudopotentials']
            self.kwargs['input_data']['input_ntyp']['Hubbard_U'] = new_data['Hubbard_U']
            self.kwargs['input_data']['input_ntyp']['starting_magnetization'] = new_data['starting_magnetization']
        else:
            pass

        atoms = self.atoms

            
        if self.step == "re-start":
            if self.sc_type == "relax":
                print("- Running relaxation step")
                calc.run_relax()
                atoms.calc = calc
                self.atoms.set_calculator(calc)
                e_relax = self.atoms.get_potential_energy()
                calc.read_results()
                self.etot[self.iteration] = e_relax
                self.opt_structures[self.iteration] = calc.atoms
            
            else:
                pass
           
            print("- Running (first) scf step")
            calc.run_scf()
            atoms.calc = calc
            atoms.set_calculator(calc)
            e_scf = atoms.get_potential_energy()
            calc.read_results()
            if self.sc_type == "scf":
                self.etot[self.iteration] = e_scf
                self.opt_structures[self.iteration] = calc.get_potential_energy()
            
            print("- Running second scf step")
            calc.run_scf2()
            if calc.scf2_necessary == True:
                atoms.set_calculator(calc)
                atoms.get_potential_energy()
            else:
                pass
            print("- Running hp step")
            calc.read_results()
            calc.run_hp()
            
            atoms_iter = calc.atoms

        elif self.step == "relax":
            print("- Restarting from (first) scf step")
            #calc.clean()
            calc.run_scf()
            atoms.calc = calc
            atoms.set_calculator(calc)
            atoms.get_potential_energy()
            calc.read_results()
                                
            self.etot[self.iteration] = calc.get_potential_energy()
            self.opt_structures[self.iteration] = calc.atoms

            
    

            print("- Running second scf step")
            calc.run_scf2()
            if calc.scf2_necessary == True:
                atoms.set_calculator(calc)
                atoms.get_potential_energy()
            else:
                pass           
 
            print("- Running hp step")
            calc.read_results()
            calc.run_hp()
            
            atoms_iter = calc.atoms

        elif self.step == "scf":
            if self.iteration ==1:
                pass
                calc.read_results()
            else:
                if self.sc_type == "relax":
                    old_directory = os.path.join(self.directory_ref,"relax_"+str(self.iteration))
                else:
                    old_directory = os.path.join(self.directory_ref,"scf_"+str(self.iteration))
                calc.read_results_fromdirectory(old_directory)
                self.etot[self.iteration] = calc.get_potential_energy()
                self.opt_structures[self.iteration] = calc.atoms

            
            print("- Restarting from hp step")
            calc.run_hp()
            atoms_iter = calc.atoms

        elif self.step == "hp":
            if self.sc_type == "relax":
                old_directory = os.path.join(self.directory_ref,"relax_"+str(self.iteration))
            else:
                old_directory = os.path.join(self.directory_ref,"scf_"+str(self.iteration))
            calc.read_results_fromdirectory(old_directory)
            self.etot[self.iteration] = calc.get_potential_energy()
            self.opt_structures[self.iteration] = calc.atoms
            atoms_iter = calc.atoms
                
           
        print("- Reading Hubbard U/V values")
        U = calc.read_Hubbard_parameters(self.iteration)
        self.UVscsd[self.iteration] =  U
        
        if self.iteration > 1:
            self.converged_U()
        else:
            pass
        print("U/V values converged: {}".format(self.converged))
        
        
        self.iteration += 1
        self.step = "re-start"
        
        new_U = update_UVscsd(U,self.Hubbard_site)
        new_data = update_Hubbard_parameters(new_U,self.kwargs['pseudopotentials'],self.Hubbard_site,atoms_iter)
        atoms_iter = new_data["atoms"]
        self.kwargs['pseudopotentials'] = new_data['pseudopotentials']
        self.kwargs['input_data']['input_ntyp']['Hubbard_U'] = new_data['Hubbard_U']
        self.kwargs['input_data']['input_ntyp']['starting_magnetization'] = new_data['starting_magnetization']
        

        
        while self.converged == False and self.iteration <= self.max_iter:
            print("-------------------------------------------")
            print("Iteration:",self.iteration,self.step)
            print("-------------------------------------------")
            calc = HpXEspresso(label=self.label_ref,
                           parallel = self.parallel,
                           parallel_hp = self.parallel_hp,
                           queue = self.queue,
                           queue_hp = self.queue_hp,
                           max_iter = self.max_iter,
                           iteration = self.iteration,
                           sc_type = self.sc_type,
                           relax_type = self.relax_type,
                           insulator = self.insulator,
                           magnetic = self.magnetic,
                           nbnd = self.nbnd,
                           parallelization_hp = self.parallelization_hp,
                           step = self.step,
                          ethr_relax = self.ethr_relax,
                           ethr_scf = self.ethr_scf,
                           ethr_scf2 = self.ethr_scf2,
                               **self.kwargs)
            

            if self.sc_type == "relax":
                print("- Running relaxation step")
                calc.run_relax()
                atoms_iter.calc = calc
                atoms_iter.set_calculator(calc)
                e_relax = atoms_iter.get_potential_energy()
                calc.read_results()
                self.etot[self.iteration] = e_relax
                self.opt_structures[self.iteration] = calc.atoms
            else:
                atoms_iter = self.atoms
                atoms_iter.calc = calc 
                atoms_iter.info['species'] = atoms_iter.get_chemical_symbols()
                for site, info in new_U.items():
                    atoms_iter.info['species'][int(site)-1] = info["new_label"]


            print("- Running (first) scf step")    
            calc.run_scf()
            atoms_iter.set_calculator(calc)
            e_scf = atoms_iter.get_potential_energy()
            calc.read_results()
            if self.sc_type == "scf":
                self.etot[self.iteration] = e_scf
            
            print("- Running second step")
            calc.run_scf2()
            atoms_iter.set_calculator(calc)
            if calc.scf2_necessary == True:
                atoms_iter.set_calculator(calc)
                atoms_iter.get_potential_energy()
            else:
                pass            

            print("- Running hp step")
            calc.read_results()
            calc.run_hp()
            
            print("- Reading Hubbard U/V values")
            U = calc.read_Hubbard_parameters(self.iteration)
            self.UVscsd[self.iteration] = U

            if self.sc_type == "relax":
                atoms_iter = calc.atoms
                
            self.converged_U()
            print("U/V values converged: {}".format(self.converged))
                
            self.iteration += 1
            self.step = "re-start"
            

            new_U = update_UVscsd(U,self.Hubbard_site)
            new_data = update_Hubbard_parameters(new_U,
                                                 self.kwargs['pseudopotentials'],
                                                 self.Hubbard_site,
                                                 atoms_iter)
            atoms_iter = new_data["atoms"]
            self.kwargs['pseudopotentials'] = new_data['pseudopotentials']
            self.kwargs['input_data']['input_ntyp']['Hubbard_U'] = new_data['Hubbard_U']
            self.kwargs['input_data']['input_ntyp']['starting_magnetization'] = new_data['starting_magnetization']
            

