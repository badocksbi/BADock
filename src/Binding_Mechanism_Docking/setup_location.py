from BioLib import *
import optparse, os, sys, subprocess

def parse_options():
    '''
    Parse arguments
    @return options object
    '''
    parser = optparse.OptionParser('setup_location.py -i PPI_NAME -j JOB_PATH -p PDB_PATH')
    parser.add_option("-i", dest="ppi_name", action="store", help="The interaction name to process [MANDATORY]", metavar="PPI_NAME")
    parser.add_option("-j", dest="job_path", action="store", help="Job path [MANDATORY]", metavar="JOB_PATH")
    parser.add_option("-p", dest="pdb_path", action="store", help="PDB path [MANDATORY]", metavar="PDB_PATH")
    (options, args) = parser.parse_args()

    if options.ppi_name == None or options.job_path == None or options.pdb_path == None:
        parser.error('Missing arguments\n')
    return options

class SetupLocation(object):
    '''
    Elucidate the location (face2face, back2back or partial) for each decoy of each potein 
    of the benckmark 5
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

        # GET AND PRINTS THE STRUCTURES FROM THE PDB FILES #
        try:
            struct_rec, struct_lig = self.__get_structures()
        except RuntimeError as e:
            sys.stdout.write(str(e))
            exit()

        # PRINTS THE DISTANCE CONSTRAINTS FILE #
        try:
            self.__setup_interface_contacts(struct_rec, struct_lig)
        except RuntimeError as e:
            sys.stdout.write(str(e))
            exit()

        sys.stdout.write('SUCCESS: Distance constraints file created!\n')
        sys.stdout.flush()

        # CHECK IF THE DIRECTED DOCKING FIND SOLUTIONS #
        self.__execute_direct_patchdock()
        direct_transform_path = os.path.join(self.job_path, 'patchdock_direct_transform.txt')
        try:
            direct_transform = self.__get_patchdock(direct_transform_path)
        except RuntimeError as e:
            sys.stdout.write(str(e))
            exit()
        if direct_transform.get_num_decoys() == 0:
            sys.stdout.write('EXIT: No decoys found in directed docking!')
            exit()

        sys.stdout.write('SUCCESS: Decoys found in direct patchdock!\n')
        sys.stdout.flush()

        # COMPUTE LOCATION FOR EACH DECOY #
        decoys_transform_path = os.path.join(self.job_path, 'patchdock_transform.txt')
        try:
            decoys_transform = self.__get_patchdock(decoys_transform_path)
        except RuntimeError as e:
            sys.stdout.write(str(e))
            exit()
        try:
            decoys_location = self.__get_location(struct_rec, struct_lig, decoys_transform)
        except RuntimeError as e:
            sys.stdout.write(str(e))
            exit()

        sys.stdout.write('SUCCESS: Location computed for each DECOY\n')
        sys.stdout.flush()

        # PRINTS THE RESULT FILE #
        self.__print_locations(decoys_location)
        sys.stdout.write('FINISHED!\n')

    def __get_structures(self):
        '''
        Parse the original PDB files in order to get the structure objects
        '''
        try:
            struct_rec = PDB.read_pdb(self.receptor_pdb, merge_chains=True)
            struct_rec.clean()
            struct_lig = PDB.read_pdb(self.ligand_pdb, merge_chains=True)
            struct_lig.clean()
        except Exception as e:
            raise RuntimeError('ERROR: Cannot parse PDB: %s\n' % str(e))
        return struct_rec, struct_lig

    def __setup_interface_contacts(self, struct_rec, struct_lig):
        '''
        Creates de distanceConstraintsFile with the interacting residues
        '''
        # Get interface residues #
        pdb_complex = Interaction(struct_rec, struct_lig)
        interacting_residues = pdb_complex.get_interacting_residues(c_type = 'CB', max_distance = 13, uniq = False)
        # Prints distanceConstraintsFile #
        dist_contr_fo = open(os.path.join(self.job_path, 'distance_constraints.txt'), 'w')
        for i in xrange(len(interacting_residues[2])):
            dist_contr_fo.write('%d %d 15\n' % (interacting_residues[0][i].get_cb().get_num(), interacting_residues[1][i].get_cb().get_num()))
        dist_contr_fo.close()

    def __get_location(self, struct_rec, struct_lig, decoys_transform):
        '''
        Returns the location of each decoy
        '''
        decoys_location = []
        for decoy in decoys_transform:
            decoy_num = decoy.get_num()
            decoy_path = os.path.join(self.job_path, 'decoys', 'decoy_%d' % decoy_num)
            decoy.print_structure(struct_rec, struct_lig, decoy_path)
            try:
                self.__execute_direct_patchdock(decoy_type='R', decoy_num=decoy_num)
                self.__execute_direct_patchdock(decoy_type='L', decoy_num=decoy_num)
            except Exception as e:
                raise e
            patchdock_R_path = os.path.join(self.job_path, 'decoys', 'decoy_%d_transformR.txt' % decoy_num)
            patchdock_L_path = os.path.join(self.job_path, 'decoys', 'decoy_%d_transformL.txt' % decoy_num)
            try:
                patchdock_R = self.__get_patchdock(patchdock_L_path)
                patchdock_L = self.__get_patchdock(patchdock_R_path)
            except RuntimeError as e:
                raise e
            num_decoys_R = patchdock_R.get_num_decoys()
            num_decoys_L = patchdock_L.get_num_decoys()
            if num_decoys_R > 0 and num_decoys_L > 0:
                location = 'E'
            elif num_decoys_R == 0 and num_decoys_L == 0:
                location = 'I'
            else:
                location = 'P'
            decoy_info = (decoy_num, location)
            decoys_location.append(decoy_info)
            sys.stdout.write('%d..' % decoy_num)
            sys.stdout.flush()
        return decoys_location

    # PATCHDOCK METHODS #

    def __execute_direct_patchdock(self, decoy_type=False, decoy_num=None):
        '''
        Execute a direct patchdock using the distance constraints
        decoy_type and decoy_num: Compute the direct patchdock between a decoy and the Ligand or Receptor
        '''
        # Patchdock execution setup #
        receptor_pdb = self.receptor_pdb
        ligand_pdb = self.ligand_pdb
        log_path = os.path.join(self.job_path, 'patchdock_direct_log.txt')
        params_path = os.path.join(self.job_path, 'patchdock_direct_params.txt')
        transform_path = os.path.join(self.job_path, 'patchdock_direct_transform.txt')
        if decoy_type == 'R':
            receptor_pdb = os.path.join(self.job_path, 'decoys', 'decoy_%d.pdb' % decoy_num)
            log_path = os.path.join(self.job_path, 'decoys', 'decoy_%d_logR.txt' % decoy_num)
            params_path = os.path.join(self.job_path, 'decoys', 'decoy_%d_paramsR.txt' % decoy_num)
            transform_path = os.path.join(self.job_path, 'decoys', 'decoy_%d_transformR.txt' % decoy_num)
        if decoy_type == 'L':
            ligand_pdb = os.path.join(self.job_path, 'decoys', 'decoy_%d.pdb' % decoy_num)
            log_path = os.path.join(self.job_path, 'decoys', 'decoy_%d_logL.txt' % decoy_num)
            params_path = os.path.join(self.job_path, 'decoys', 'decoy_%d_paramsL.txt' % decoy_num)
            transform_path = os.path.join(self.job_path, 'decoys', 'decoy_%d_transformL.txt' % decoy_num)
        params_fo = open(params_path, 'w')
        params_fo.write(self.__create_patchdock_params(receptor_pdb, ligand_pdb, log_path))
        params_fo.close()
        # Execute patchdock #
        patchdock_cmd = '/soft/bio/structure/PatchDock/patch_dock.Linux %s %s' % (params_path, transform_path)
        try:
            subprocess.check_output(patchdock_cmd, stderr=subprocess.STDOUT, shell=True)
        except Exception as e:
            raise RuntimeError('ERROR: Error executing PatchDock: %s\n' % str(e))
        os.remove(params_path)

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

    def __create_patchdock_params(self, receptor_pdb, ligand_pdb, log_file):
        '''
        Returns the patchdock params string
        '''
        params =  'receptorPdb %s\nligandPdb %s\n' % (receptor_pdb, ligand_pdb)
        params += 'distanceConstraintsFile %s\n' % os.path.join(self.job_path, 'distance_constraints.txt')
        params += 'protLib /soft/bio/structure/PatchDock/chem.lib\n'
        params += 'log-file %s\nlog-level 2\n' % log_file
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

    def __print_locations(self, decoys_location):
        '''
        Prints the location for each decoy
        '''
        result_fo = open(os.path.join(self.job_path, 'decoys_location.txt'), 'w')
        decoys_location.sort()
        for decoy in decoys_location:
            result_fo.write('%d\t%s\n' % decoy)
        result_fo.close()

#--------------------------#
#          MAIN            #
#--------------------------#

if __name__ == "__main__":

    options = parse_options()
   
    setup_location = SetupLocation(options.ppi_name, options.job_path, options.pdb_path) 
    setup_location.run()