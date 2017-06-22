from BioLib import *
import optparse, os, sys, subprocess

def parse_options():
    '''
    Parse arguments
    @return options object
    '''
    parser = optparse.OptionParser('setup_decoys.py -i PPI_NAME -j JOB_PATH -p PDB_PATH')
    parser.add_option("-i", dest="ppi_name", action="store", help="The interaction name to process [MANDATORY]", metavar="PPI_NAME")
    parser.add_option("-j", dest="job_path", action="store", help="Job path [MANDATORY]", metavar="JOB_PATH")
    parser.add_option("-p", dest="pdb_path", action="store", help="PDB path [MANDATORY]", metavar="PDB_PATH")
    (options, args) = parser.parse_args()

    if options.ppi_name == None or options.job_path == None or options.pdb_path == None:
        parser.error('Missing arguments\n')
    return options

class SetupDecoys(object):
    '''
    Create a set of decoys of known ppi (benchmark5) using rigit body docking (PatchDock) and
    refinement (FiberDock) in order to further study the different type of interfaces (face2face, 
    back2back and face2back)
    '''
    def __init__(self, ppi, job_path, pdb_path):
        '''
        Contructor
        '''
        self.ppi = ppi
        self.job_path = job_path
        self.pdb_path = pdb_path

        self.receptor_pdb = os.path.join(self.pdb_path, self.ppi+'_r_b.pdb')
        self.ligand_pdb = os.path.join(self.pdb_path, self.ppi+'_l_b.pdb')

    def run(self):
        '''
        Start the execution
        '''
        sys.stdout.write('--> START %s\n' % self.ppi)
        sys.stdout.flush()

        # DOCK UNBOUND RECEPTOR AND LIGAND WITH PATCHDOCK TO GET THE DECOYS #
        try:
            self.__execute_patchdock()
        except RuntimeError as e:
            sys.stdout.write(str(e))
            exit()
        sys.stdout.write('SUCCESS: Decoys Computed!\n')
        sys.stdout.flush()

        # REFINE AND PRINT SRUCTURES WITH FIBERDOCK #
        try:
            self.__execute_fiberdock()
        except RuntimeError as e:
            sys.stdout.write(str(e))
            exit()
        sys.stdout.write('SUCCESS: Decoys Refined and Printed!\n')
        sys.stdout.flush()

        sys.stdout.write('FINISHED!\n')

    def __execute_patchdock(self):
        '''
        Execute patchdock
        '''
        # Patchdock execution setup #
        patchdock_params_path = os.path.join(self.job_path, 'patchdock_params.txt')
        transform_path = os.path.join(self.job_path, 'patchdock_transform.txt')
        patchdock_params_fo = open(patchdock_params_path, 'w')
        patchdock_params_fo.write(self.__create_patchdock_params())
        patchdock_params_fo.close()
        # Execute patchdock #
        patchdock_cmd = '/soft/bio/structure/PatchDock/patch_dock.Linux %s %s' % (patchdock_params_path, transform_path)
        try:
            subprocess.check_output(patchdock_cmd, stderr=subprocess.STDOUT, shell=True)
        except Exception as e:
            raise RuntimeError('ERROR: Error executing PatchDock: %s\n' % str(e))
        # Rewrite patchdock transform output and include native as solution 0 #
        try:
            patchdock_decoys = self.__get_patchdock(transform_path)
        except RuntimeError as e:
            raise e
        transform_clean_path = os.path.join(self.job_path, 'patchdock_transform.clean.txt')
        transform_clean_fo = open(transform_clean_path, 'w')
        for patchdock_decoy in patchdock_decoys:
            transform_clean_fo.write('%d %f %f %f %f %f %f\n' % (patchdock_decoy.get_num(), 
                                                                 patchdock_decoy.alpha,
                                                                 patchdock_decoy.beta,
                                                                 patchdock_decoy.gamma,
                                                                 patchdock_decoy.x,
                                                                 patchdock_decoy.y,
                                                                 patchdock_decoy.z))
        transform_clean_fo.write('%d 0 0 0 0 0 0\n' % (patchdock_decoys.get_num_decoys()+1))
        transform_clean_fo.close()
        os.remove(patchdock_params_path)

    def __get_patchdock(self, patchdock_path):
        '''
        PARSE PATCHDOCK transformation file
        '''
        if not os.path.isfile(patchdock_path):
            raise RuntimeError('ERROR: No PATCHDOCK TRASNFORMATION file found\n')
        try:
            patchdock_results = PATCHDock(patchdock_path)
        except Exception as e:
            raise RuntimeError('ERROR: Cannot parse PATCHDOCK TRASNFORMATION file: %s\n' % str(e))
        return patchdock_results

    def __create_patchdock_params(self):
        '''
        Returns the patchdock params string
        '''
        params =  'receptorPdb %s\nligandPdb %s\n' % (self.receptor_pdb, self.ligand_pdb)
        params += 'protLib /soft/bio/structure/PatchDock/chem.lib\n'
        params += 'log-file %s\nlog-level 2\n' % os.path.join(self.job_path, 'patchdock_log.txt')
        params += 'receptorSeg 10.0 20.0 1.5 1 0 1 0\n'
        params += 'ligandSeg 10.0 20.0 1.5 1 0 1 0\n'
        params += 'scoreParams 0.3 -5.0 0.5 0.0 0.0 1500 -8 -4 0 1 0\n'
        params += 'desolvationParams 500.0 1.0\n'
        params += 'clusterParams 0.1 4 2.0 9.0\n'
        params += 'baseParams 4.0 13.0 2\n'
        params += 'matchingParams 1.5 1.5 0.4 0.5 0.9\n'
        params += 'matchAlgorithm 1\n'
        params += 'receptorGrid 0.5 6.0 6.0\n'
        params += 'ligandGrid 0.5 6.0 6.0\n'
        params += 'receptorMs 10.0 1.8\n'
        params += 'ligandMs 10.0 1.8\n'
        return params

    def __execute_fiberdock(self):
        '''
        Execute FiberDock
        '''
        # Create PDBs with Hydrogens #
        add_hydrogens_cmd = '/soft/bio/structure/FiberDock1.1/addHydrogens.pl %s %s' % (self.receptor_pdb, self.ligand_pdb)
        try:
            #subprocess.check_output(add_hydrogens_cmd, stderr=subprocess.STDOUT, shell=True)
            pass
        except Exception as e:
            raise RuntimeError('ERROR: Error executing addHydrogens: %s\n' % str(e))
        # Calculate normal modes for backbone refinement using only CA #
        only_CA_cmd = 'cat %s | grep ATOM | grep " CA " > %s.CA'
        try:
            subprocess.check_output(only_CA_cmd % (self.receptor_pdb, self.receptor_pdb), stderr=subprocess.STDOUT, shell=True)
            subprocess.check_output(only_CA_cmd % (self.ligand_pdb, self.ligand_pdb), stderr=subprocess.STDOUT, shell=True)
        except Exception as e:
            raise RuntimeError('ERROR: Error creating CA PDBs: %s\n' % str(e))
        nma_cmd = '/soft/bio/structure/FiberDock1.1/nma %s.CA %s.nma 85 3 10'
        try:
            subprocess.check_output(nma_cmd % (self.receptor_pdb, os.path.join(self.job_path, self.ppi+'_r_b')), stderr=subprocess.STDOUT, shell=True)
            subprocess.check_output(nma_cmd % (self.ligand_pdb, os.path.join(self.job_path, self.ppi+'_l_b')), stderr=subprocess.STDOUT, shell=True)
        except Exception as e:
            raise RuntimeError('ERROR: Error executing nma: %s\n' % str(e))
        # Execute FiberDock #
        fiberdock_params_path = os.path.join(self.job_path, 'fiberdock_params.txt')
        fiberdock_params_fo = open(fiberdock_params_path, 'w')
        fiberdock_params_fo.write(self.__create_fiberdock_params())
        fiberdock_params_fo.close()
        fiberdock_exe_cmd = '/soft/bio/structure/FiberDock1.1/FiberDock %s' % fiberdock_params_path
        try:
            subprocess.check_output(fiberdock_exe_cmd, stderr=subprocess.STDOUT, shell=True)
        except Exception as e:
            raise RuntimeError('ERROR: Error executing FiberDock: %s\n' % str(e))
        os.remove(fiberdock_params_path)

    def __create_fiberdock_params(self):
        '''
        Returns the fiberdock params string
        '''
        ### I/O
        params =  'receptorPDBFileName %s.HB\n' %  self.receptor_pdb
        params += 'ligandPDBFileName %s.HB\n' % self.ligand_pdb
        ## reference for rmsd calculations
        #params += 'templateLigandPDBFileName %s.HB\n' % self.ligand_pdb
        # transformations for refinement
        params += 'transFileName %s\n' % os.path.join(self.job_path, 'patchdock_transform.clean.txt')
        ## libraries files
        params += 'rotamerLibFile  /soft/bio/structure/FiberDock1.1/lib/bbdep02.May.sortlib\n'
        params += 'protLib /soft/bio/structure/FiberDock1.1/lib/chem.lib\n'
        ## output file
        params += 'energiesOutFileName %s\n' % os.path.join(self.job_path, 'fiberdock_energies')
        ### Output Options
        # to output refined complexes
        params += 'printRefinedComplexes 0\n'
        ## 1 - only energy caclulaltion is performed without refinement (works only for FiberDock pre)
        params += 'onlyEnergyCalculation 0\n'
        ### side-chain optimization
        # 1 - only clashing residues are flexible, 0 - all residues are flexible
        params += 'receptorOnlyClashesMovable 1\n'
        params += 'ligandOnlyClashesMovable 1\n'
        # 0 - small rotamer set, 1 - extended rotamer set
        params += 'extraRotamers 0\n'
        params += 'lpMethod 2\n'
        params += 'exclusionsFile /soft/bio/structure/FiberDock1.1/lib/exclusions\n'
        ### rigid-body optimization (RBO)
        # num of MC cycles for rigid-body minimization (if 0 - no RBO)
        params += 'rigidBodyMinimizationCycles 10\n'
        params += 'charmmParamFile /soft/bio/structure/FiberDock1.1/lib/par_all27_prot_na.prm\n'
        params += 'definitionsFile /soft/bio/structure/FiberDock1.1/lib/top_all27_prot_na.rtf\n'
        ### weights for energy score (Default)
        params += 'attrVdWWeight   1.5\n'
        params += 'repVdWWeight    0.8\n'
        params += 'ACEWeight       1.6\n'
        params += 'attrElWeight    0.21\n'
        params += 'repElWeight     0.21\n'
        params += 'l_attrElWeight  0.0\n'
        params += 'l_repElWeight   0.69\n'
        params += 'HBWeight        1.2\n'
        params += 'pipiWeight      1.0\n'
        params += 'catpiWeight     0.7\n'
        params += 'aliphWeight     2.5\n'
        params += 'insidenessWeight        0.7\n'
        params += 'confProbWeight  0.0\n'
        params += 'radiiScaling 0.80\n'
        ### Backbone refinement parameters
        params += 'performBackboneRef 1\n'
        params += 'Cycles 10\n'
        params += 'RecAnalyzedFile %s.CA\n' % self.receptor_pdb
        params += 'ReceptorNMs %s_r_b.nma\n' % os.path.join(self.job_path, self.ppi)
        params += 'Lambda 0.05\n'
        params += 'LigAnalyzedFile %s.CA\n' % self.ligand_pdb
        params += 'LigandNMs %s_l_b.nma\n' % os.path.join(self.job_path, self.ppi)
        ## LOG-FILE params
        params += 'log-file %s\n' % os.path.join(self.job_path, 'fiberdock_log.txt')
        params += 'log-level 2\n'
        return params

#--------------------------#
#          MAIN            #
#--------------------------#

if __name__ == "__main__":

    options = parse_options()
   
    setup_decoys = SetupDecoys(options.ppi_name, options.job_path, options.pdb_path) 
    setup_decoys.run()