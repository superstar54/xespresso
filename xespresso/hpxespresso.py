from ase import Atoms
from ase import io
from ase.build import bulk, sort
from ase.io import read
from xespresso import Espresso
from itertools import chain
import numpy as np
import copy
import os, shutil
import sys

class HpXEspresso(Espresso):
    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='xespresso', atoms=None, package = 'pw', parallel =  '',
                 queue = None,queue_hp = False, parallel_hp = '', iteration=1,
                 sc_type ="relax", relax_type = "vc-relax",
                 insulator = None, magnetic = None, nbnd =  None,
                 parallelization_hp = False,ethr_relax = None, ethr_scf = None, ethr_scf2 = None,
                 **kwargs):
        """
        This class allows to perform the chain of relax-scf-hp calculations necessary to compute U/V values
        via DFPT "Hubbard parameters from density-functional perturbation theory", 
        Phys. Rev. B 98, 085127 (2018); arXiv:1805.01805 using ase and XEspresso, 
        of which XpXEspresso is a subclass.
        
        input_data, pseudopotentials, kspacing, kpts, koffset
            Please have a look at Espresso module in ASE
            
        queue: dict
            A dictionary with parameters for job submission of the pw.x calculation, e.g.
             queue = {'nodes': 4, 'ntasks-per-node': 20, 
                      'account': 'xxx', 'partition': 'normal', 
                      'time': '23:59:00'}
                      
        queue_hp: dict
            A dictionary with parameters for job submission of the hp.x calculation, e.g.
             queue = {'nodes': 4, 'ntasks-per-node': 20, 
                      'account': 'xxx', 'partition': 'normal', 
                      'time': '23:59:00'}
        package: str
            Choose the quantum espresso pacakge: pw, dos, projwfc, band, pp, ph, ..
            For NEB calculation, please use neb.NEBEspresso module.
            Calculaiton use phonon is not implemented yet.
            
        parallel: str
            A str which control the parallelization parameters for the pw.x runs: -nimage, -npools, 
            -nband, -ntg, -ndiag or -northo (shorthands, respectively: 
            -ni, -nk, -nb, -nt, -nd).
        
        parallel_hp: str
            A str which control the parallelization parameters for the hp.x runs: -nimage, -npools, 
            -nband, -ntg, -ndiag or -northo (shorthands, respectively: 
            -ni, -nk, -nb, -nt, -nd).
            
        max_iter: int
            An int indicating the interation step in case a self-consistent (site-dependent) calculation
            of the U/V parameters is performed (Chech the Uscsd module).
        
        sc_type: str 
            A str indicating if the procedure for the self-consistent (site-dependent) calculation
            of the U/V parameters is performed by relaxing or not the structure. Check also the Uscsd class.
            Possible values are "relax" and "scf".
            
        relax_type: str
            Possible values are "vc-relax" and "relax". Indicate thetype of relaxation run incase sc_type = "relax".
            Check also the Uscsd class.
        
        insulator: bool
            A bool indicating if the material is an insulator (True) or a metal (False).
            It has to be specified when sc_type = "scf" and iteration = 1, otherwhise one could not determine 
            the correct procedure to be perfomed.
            
        magnetic: int
            An int indicating if the material is non-magneticc (1) or  magnetic (2)
            It has to be specified when sc_type = "scf" and iteration = 1, otherwhise one could not determine 
            the correct procedure to be perfomed.
            
        nbnd: int
            An int indicating the number of bands to take into account in the calculations.
            It has to be specified when sc_type = "scf" and iteration = 1, otherwhise one could not determine 
            the correct procedure to be perfomed.
            
        parallelization_hp: bool
            A bool indicating if  parallelization of the hp run over atoms or qpoints needs to be performed.
            parallelization_hp = True
            
        ethr_relax: float
            A float corresponfing to the threshold for the scf convergence in the relaxation step.  
            If not specified it is set to 1.0E-8
        
        ethr_scf: float
            A float corresponfing to the threshold for the scf convergence in the (first) scf step.
            If not specified it is set to 1.0E-10
            
        ethr_scf2: float
            A float corresponfing to the threshold for the scf convergence in the second scf step.
            If not specified it is set to 1.0E-12
            
        !!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!
        Do not forget to provide in the input_data dictionry and two elements:
            "lda_plus_U" : True
            "input_ntyp": { "starting_magnetization" : {Hubbard_site: XX},
                            "Hubbard_U" : {Hubbard_site: XX}}
            
        TO DO:
        - check inputs: sc_type, relax_typ, insulator, magnetic, ndnb
        - implement parallelization_hp = True
        - implement possibility of restart of an interrupted relax run or hp run.
        - check copy parameters.in file in the case of DFT+U+V
            
        """
        super().__init__(restart=None, ignore_bad_restart_file=False,
                 label='xespresso', atoms=None, package = 'pw', parallel = '',
                 queue = None, 
                 **kwargs)
        
        self.calc = None
        
        self.relax_directory = None
        self.relax_parameters = None
        self.relax_results = None
        
        self.scf_directory = None
        self.scf_parameters = None
        self.scf_results = None
        
        self.hp_directory = None
        if 'parameters_hp' not in kwargs:
            kwargs['parameters_hp'] = {}
        self.hp_parameters = kwargs['parameters_hp']
        self.hp_results = None
        
        self.directory = os.path.split(label)[0]
        self.prefix = os.path.split(label)[1]
        self.directory_ref = self.directory
        

        self.queue_hp = queue_hp
        self.parallel_hp = parallel_hp
        self.queue = queue
        self.parallel = parallel
            

        self.iteration = iteration
        
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
        
        
        self.ethr_relax = ethr_relax
        self.ethr_scf = ethr_scf
        self.ethr_scf2 = ethr_scf2
            
    def read_results_fromdirectory(self,directory):
        """
        Extension of the read_results metho of the XEspresso class, so that results
        can be read from a directory different from self.directory
        
        directory: str
            A str indicating the the path of the directory from which results have to be read.
        """
        from xespresso.xio import get_atomic_species
        pwo = self.prefix + '.pwo'
        
        if directory == None:
            directory = self.directory
        
        pwos = [file for file in os.listdir(directory) if pwo in file]
        output = None
        for pwo in pwos:
            atomic_species = None
            pwo = os.path.join(directory, pwo)
            # print('read results from: %s'%pwo)
            # output = next(read_espresso_out(pwo))
            atomic_species = get_atomic_species(pwo)
            try:
                output = io.read(pwo)
                atomic_species = get_atomic_species(pwo)
                if atomic_species: output.info['species'] = atomic_species
            except Exception as e:
                print(pwo, e)
                # print('\nread %s failed\n'%pwo)
        if output:
            self.calc = output.calc
            self.results = output.calc.results
            self.results['atoms'] = output
            self.efermi = self.get_fermi_level()
            self.nspins = self.get_number_of_spins()
            self.atoms = output
        else:
            print('\nRead result failed\n') 
        
        
    def is_insulator(self,fermi=None):
        """
        Checking if the system is a insulators
        One can provide the fermi energy as input, otherwise it will be extracted from the pw.pwo file.
        WARNING: The Fermi energy cannot beread from a calculation run using "tot_magnetization"
        
        fermi: float
            A float corresponding to the Fermi energy
        
        Results:
        - set with Bool (True=insulator, False=metal), bandgap value (eV), VBM and CBm values (eV)
        """

        if fermi == None:
            fermi = self.calc.get_fermi_level()
            if fermi == None:
                raise Exception("I cannot get the Fermi energy from the pw.pwo file"
                                "Provide a Fermi energy value.")


        # reorganize the bands, rather than per kpoint, per energy level
        bands = []
        nbands = len(self.calc.get_eigenvalues(kpt=0,spin=0))
        for i in range(nbands):
            bands.append([])
            for k in range(len(self.calc.get_k_point_weights())):
                if self.calc.get_number_of_spins() == 1:#calc.nspins == 1:
                    bands[i].append(self.calc.get_eigenvalues(kpt=k)[i])
                else:
                    bands[i].append(self.calc.get_eigenvalues(kpt=k,spin=0)[i])
                    bands[i].append(self.calc.get_eigenvalues(kpt=k,spin=1)[i])
        bands = np.array(bands)
        max_mins = [(max(i), min(i)) for i in bands]

        # one band is crossed by the fermi energy
        if any(i[1] < fermi and fermi < i[0] for i in max_mins):
            return(False, None, None, None)

        # case of semimetals, fermi energy at the crossing of two bands
        # this will only work if the dirac point is computed!
        elif (any(i[0] == fermi for i in max_mins) and
                  any(i[1] == fermi for i in max_mins)):
            return False, 0., fermi, fermi
        # insulating case
        else:
            # take the max of the band maxima below the fermi energy
            homo = max([i[0] for i in max_mins if i[0] < fermi])
            # take the min of the band minima above the fermi energy
            lumo = min([i[1] for i in max_mins if i[1] > fermi])
            gap = lumo - homo
            if gap <= 0.:
                raise Exception("Something wrong has been implemented. "
                                "Revise the code!")
            return True, gap, homo, lumo
        
    def get_magnetizations(self, getall=False):
        """
        Read the total and absolute magnetizations from the pw.x out file
        """

        start = '!    total energy'
        end = '     convergence has'

        mags = dict()
        mags['total magnetization'] = []
        mags['absolute magnetization'] = []
        
        filename = os.path.join(self.directory, self.prefix + ".pwo")
        with open(filename, 'r') as fobj:
            for line in fobj:
                if start in line:
                    while end not in line:
                        line = next(fobj)
                        for mag in mags.keys():
                            if mag in line:
                                mags[mag].append(float(line.split()[3]))
        if getall:
            return mags
        else:
            return [mags[k][-1] if mags[k] else 0.0 for k in
                    sorted(mags.keys())]
        
    def read_Hubbard_parameters(self,iteration):
        """
        Read the U/V values from the Hubbard_parameters.dat file or from the paramets.in file for DFT+U+V
        """
        
        
        if 'lda_plus_u_kind' not in self.parameters or self.parameters['lda_plus_u_kind'] != 1:
            start = 'site n.'
            hp_U = dict()
            directory = os.path.join(self.directory_ref,"hp_"+str(iteration))
            filename = os.path.join(directory, self.prefix + ".Hubbard_parameters.dat")
            
            print(filename)
            try:
                file = open(filename, 'r')
            except IOError:
                print('There was an error opening the Hubbard.dat file!')
                sys.exit('There was an error opening the Hubbard.dat file!')

            
            with open(filename, 'r') as fobj:
                lines = fobj.readlines()
                nlines = len(lines)
                start = [n  for n, line in enumerate(lines) if 'site n.' in line][0]
                end = [n  for n, line in enumerate(lines) if "=--------------" in line][1] -1
                for line in lines[start+1:end]:
                    l = line.split()
                    hp_U[l[0]] = {}
                    hp_U[l[0]]['label'] = l[2]
                    hp_U[l[0]]['spin'] = l[3]
                    hp_U[l[0]]['type'] = l[4]
                    hp_U[l[0]]['new_label'] = l[5]
                    hp_U[l[0]]['U'] = l[6]

        #elif self.scf_parameters['lda_plus_u_kind'] == 2:
        #    hp_U = None
        #    print("Reading the parametes.in file for DFT+U+V is not implemented yet")
        else:
            hp_U = None
            print("lda_plus_u_kind=1 is not implemented")

        return hp_U
        
    def run_relax(self):
        """
        Running the relaxation step
        """
        
        if self.sc_type == "scf":
            sys.exit("You set sc_type to False. I will not relax the structure")
        
        self.directory = os.path.join(self.directory_ref,"relax_"+str(self.iteration))
        self.set_queue(package = "pw", parallel = self.parallel, queue = self.queue)
        
        if self.relax_parameters != None: # This works if the execution of the script is not interrupted
            self.parameters = copy.deepcopy(self.relax_parameters)
        
        self.parameters["calculation"] = self.relax_type
        self.parameters['verbosity'] = 'high'
        if self.ethr_relax == None:
            self.parameters["conv_thr"] = 1.0E-8
        else:
            self.parameters["conv_thr"] = self.ethr_relax
            
        self.parameters['input_data']['restart_mode'] = 'from_scratch'
        
        if self.iteration > 1:
            if 'lda_plus_u_kind' in self.parameters['input_data'] and self.parameters['input_data']['lda_plus_u_kind'] == 2:
                old_parametersin =  os.path.join(self.directory_ref,"hp_"+str(self.iteration -1),'parameters.in')
                shutil.copy(old_parametersin,
                            self.directory ,symlinks=False, ignore=None)
                self.parameters['input_data']['Hubbard_parameters'] = 'file'

    def run_scf(self):
        """
        Running the (first) scf step 
        """
        #self.set_queue(package = "pw", parallel = self.parallel, queue = self.queue)
        if self.sc_type == "relax":
            #self.read_results()
            if self.iteration == 1:
                self.read_results()
            else:
                previous_run = os.path.join(self.directory_ref,"relax_"+str(self.iteration))
                self.read_results_fromdirectory(directory=previous_run)
                self.relax_directory = previous_run
                
            if not self.relax_directory:
                self.relax_directory = self.directory
            if not self.relax_parameters:
                self.relax_parameters = copy.deepcopy(self.parameters)
            if not self.relax_results:
                self.relax_results = copy.deepcopy(self.results)
                self.efermi = self.get_fermi_level()
            # read new atoms from results and set magmoms
            self.atoms = self.results['atoms']
            self.parameters = copy.deepcopy(self.relax_parameters)
            self.parameters['input_data']['restart_mode'] = 'restart'
            
            
            # create working directory
            self.save_directory_old = os.path.join(self.relax_directory, '%s.save'%self.prefix)
            self.directory = os.path.join(self.directory_ref,"scf_"+str(self.iteration))
            self.save_directory = os.path.join(self.directory, '%s.save'%self.prefix)
            self.label = os.path.join(self.directory, self.prefix)
            self.scf_directory = self.directory
            
            if self.iteration > 1:
                if 'lda_plus_u_kind' in self.parameters['input_data'] and self.parameters['input_data']['lda_plus_u_kind'] == 2:
                    old_parametersin =  os.path.join(self.directory_ref,"hp_"+str(self.iteration -1),'parameters.in')
                    shutil.copy(old_parametersin,
                                self.directory ,symlinks=False, ignore=None)
                    self.parameters['input_data']['Hubbard_parameters'] = 'file'
    


            #if path already exists, remove it before copying with copytree()
            if os.path.exists(self.save_directory):
                shutil.rmtree(self.save_directory)
            
            shutil.copytree(self.save_directory_old,
                            self.save_directory,symlinks=False, ignore=None)
            #self.calculate()
            insulator = self.is_insulator()[0]
            magnetic = self.get_number_of_spins()
            nbnd = len(self.calc.get_eigenvalues(kpt=0))
        else:                                  
            self.directory = os.path.join(self.directory_ref,"scf_"+str(self.iteration))
            self.scf_directory = self.directory
                                                 
            if self.iteration == 1:
                if self.insulator == None:
                    sys.exit("You need to provide a value (True/False) for  'insulator'")
                else:
                    insulator = self.insulator
                if self.magnetic == None:
                    sys.exit("You need to provide a value (True/False) for  'magnetic'")
                else:
                    magnetic = self.magnetic
                if self.nbnd == None:
                    sys.exit("You need to provide a value (True/False) for  'nbnd'")
                else:
                    nbnd = self.nbnd
                    
            else:
                previous_run = os.path.join(self.directory_ref,"scf_"+str(self.iteration-1))
                self.read_results_fromdirectory(directory=previous_run)
                insulator = self.is_insulator()[0]
                magnetic = self.get_number_of_spins()
                nbnd = len(self.calc.get_eigenvalues(kpt=0))
                
       
                
        self.parameters['verbosity'] = 'high'
        self.parameters['input_data']['verbosity'] = 'high'
        self.parameters['calculation'] = "scf"
        
        if insulator == True and magnetic==1:
            self.parameters['input_data']["occupations"] = "fixed"
            #nbnd = len(self.calc.get_eigenvalues(kpt=0))
            self.parameters['input_data']["nbnd"] = nbnd
            if  "degauss" in self.parameters['input_data'].keys():
                del self.parameters['input_data']["degauss"]

            if  "starting_magnetization" in self.parameters['input_data']['input_ntyp'].keys():
                del self.parameters['input_data']['input_ntyp']["starting_magnetization"]
            
            if self.ethr_scf == None:
                self.parameters["conv_thr"] = 1.0E-10
            else:
                self.parameters["conv_thr"] = self.ethr_scf


        elif insulator == False or (insulator == True and magnetic==2):
            self.parameters['input_data']["occupations"] = "smearing"
            self.parameters['input_data']["degauss"] = 0.01
            
            if self.ethr_scf == None:
                self.parameters["conv_thr"] = 1.0E-10
            else:
                self.parameters["conv_thr"] = self.ethr_scf
          
        else:
            if self.sc_type =="relax":
                sys.exit("Case not implemented: is an insulator? {} and is spin polarized? {}".format(self.is_insulator()[0],self.get_spin_polarized()) )
            else:
                sys.exit("Case not implemented: is an insulator? {} and is spin polarized? {}".format(self.insulator,self.magnetic))

        self.set_queue(package = "pw", parallel = self.parallel, queue = self.queue)
    
    def run_scf2(self):
        """
        Running the second scf step (if necessary)
        """
        self.read_results()       
        if self.get_number_of_spins() == 2 and self.is_insulator()[0] == True:

                # save relax parameters
                if not self.scf_directory:
                    self.scf_directory = self.directory
                if not self.scf_parameters:
                    self.scf_parameters = copy.deepcopy(self.parameters)
                if not self.scf_results:
                    self.scf_results = copy.deepcopy(self.results)
                    self.efermi = self.get_fermi_level()

                # read new atoms from results and set scf specific parameters 
                self.atoms = self.results['atoms']
                self.parameters = copy.deepcopy(self.scf_parameters)
                self.parameters['calculation'] = "scf"
                self.parameters['input_data']['verbosity'] = 'high'
                self.parameters['input_data']['restart_mode'] = 'restart'
                self.parameters['input_data']["occupations"] = "fixed"
                nbnd = len(self.calc.get_eigenvalues(kpt=0))
                self.parameters['input_data']["nbnd"] = nbnd
                self.parameters['input_data']["tot_magnetization"] = abs(self.get_magnetizations(True)['total magnetization'][-1])

                if  "degauss" in self.parameters['input_data'].keys():
                    del self.parameters['input_data']["degauss"]
                    
                if  "starting_magnetization" in self.parameters['input_data']['input_ntyp'].keys():
                    del self.parameters['input_data']['input_ntyp']["starting_magnetization"]
                    
                self.parameters['input_data']["startingpot"] = "file"
                self.parameters['input_data']["startingwfc"] = "file"

                if self.ethr_scf2 == None:
                    self.parameters["conv_thr"] = 1.0E-12
                else:
                    self.parameters["conv_thr"] = self.ethr_scf2
                


                # create working directory
                self.save_directory_old = os.path.join(self.scf_directory, '%s.save'%self.prefix)
                self.directory = os.path.join(self.directory_ref,"scf2_"+str(self.iteration))
                self.save_directory = os.path.join(self.directory, '%s.save'%self.prefix)
                self.label = os.path.join(self.directory, self.prefix)
                self.scf_directory = self.directory
                
                #if path already exists, remove it before copying with copytree()
                if os.path.exists(self.save_directory):
                    shutil.rmtree(self.save_directory)

                shutil.copytree(self.save_directory_old,
                                self.save_directory ,symlinks=False, ignore=None)
                
                self.scf2_necessary = True
        else:
            self.scf2_necessary = False
            print("No need to perform a second scf")
            pass
   
        self.set_queue(package = "pw", parallel = self.parallel, queue = self.queue)
            
    def run_hp(self): 
        """
        Running the hp step
        """
        # create working directory
        if self.scf_directory == None:
            if os.path.exists(os.path.join(self.directory_ref,"scf2_"+str(self.iteration))):
                self.scf_directory = os.path.join(self.directory_ref,"scf2_"+str(self.iteration))
            else:
                self.scf_directory = os.path.join(self.directory_ref,"scf_"+str(self.iteration))
          
        if self.iteration == 1:
            self.read_results()  
        else:
            self.read_results_fromdirectory(directory=self.scf_directory)


            
        if self.parallelization_hp == False:
            if not self.hp_parameters:
                self.hp_parameters = copy.deepcopy(self.parameters_hp)
            if not self.hp_results:
                self.hp_results = copy.deepcopy(self.results)

            
            self.save_directory_old = os.path.join(self.scf_directory, '%s.save'%self.prefix)
            self.directory = os.path.join(self.directory_ref,"hp_"+str(self.iteration))
            #self.directory = os.path.join(self.directory_ref,"hp_",str(iteration))
            self.save_directory = os.path.join(self.directory, '%s.save'%self.prefix)
            self.label = os.path.join(self.directory, self.prefix)
            self.hp_directory = self.directory

            #if path already exists, remove it before copying with copytree()
            if os.path.exists(self.save_directory):
                shutil.rmtree(self.save_directory)

            shutil.copytree(self.save_directory_old,
                            self.save_directory ,symlinks=False, ignore=None)

            self.post(package = 'hp',queue = self.queue_hp, parallel = self.parallel_hp, **self.hp_parameters)
        else:
            sys.exit("Parallelization of the hp.x calculation over the atoms or q-points not implemented.")
            
            

