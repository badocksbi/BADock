from BioLib import *
import optparse, os, sys, subprocess

def parse_options():
    '''
    Parse arguments
    @return options object
    '''
    parser = optparse.OptionParser('setup_score.py -i PPI_NAME -j JOB_PATH -p PDB_PATH')
    parser.add_option("-i", dest="ppi_name", action="store", help="The interaction name to process [MANDATORY]", metavar="PPI_NAME")
    parser.add_option("-j", dest="job_path", action="store", help="Job path [MANDATORY]", metavar="JOB_PATH")
    parser.add_option("-p", dest="pdb_path", action="store", help="PDB path [MANDATORY]", metavar="PDB_PATH")
    (options, args) = parser.parse_args()

    if options.ppi_name == None or options.job_path == None or options.pdb_path == None:
        parser.error('Missing arguments\n')
    return options

class SetupScore(object):
    '''
    Score each PPI
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

        # GET THE STRUCTURES FROM THE PDB FILES #
        try:
            struct_rec, struct_lig = self.__get_structures()
        except RuntimeError as e:
            sys.stdout.write('[\033[91mFAIL\033[0m]' % str(e))
            exit()

        # COMPUTE ILOOPS AND SPLITPOTENTIALS SCORES FOR EACH DECOY #     
        decoys_transform_path = os.path.join(self.job_path, 'patchdock_transform.txt')
        try:
            decoys_transform = self.__get_patchdock(decoys_transform_path)
        except RuntimeError as e:
            sys.stdout.write(str(e))
            exit()
        decoy_scores = self.__calculate_decoy_scores(struct_rec, struct_lig, decoys_transform)
        sys.stdout.write('[\033[92mSUCCESS\033[0m] iLoops and SplitPotentials computed\n')
        sys.stdout.flush()

        # PRINT SCORES #
        self.__print_energies(decoy_scores)

        sys.stdout.write('[\033[92mFINISHED\033[0m]\n')

    def __get_structures(self):
        '''
        Parse the original PDB files in order to get the structure objects
        '''
        try:
            struct_rec = PDB.read_pdb(self.receptor_pdb, merge_chains=True)
            struct_rec.clean()
            struct_rec.set_dssp()
            struct_lig = PDB.read_pdb(self.ligand_pdb, merge_chains=True)
            struct_lig.clean()
            struct_lig.set_dssp()
        except Exception as e:
            raise RuntimeError('ERROR: Cannot parse PDB: %s\n' % str(e))
        return struct_rec, struct_lig

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

    def __calculate_decoy_scores(self, struct_rec, struct_lig, decoys_transform):
        '''
        Calculates for each decoy the split potentials and iLoops scores.
        '''
        decoy_scores = []
        for decoy in decoys_transform:
            decoy_structure = decoy.get_structure(struct_lig)
            decoy_interaction = Interaction(struct_rec, decoy_structure)
            try:
                splitpotentials_energies = self.__get_splitpotentials(decoy_interaction)
            except RuntimeError as e:
                sys.stdout.write(str(e))
                sys.stdout.flush()
                continue
            rmsd = struct_lig.get_RMSD(decoy_structure)
            decoy_scores.append((decoy.get_num(), rmsd, splitpotentials_energies))
        return decoy_scores

    def __get_splitpotentials(self, decoy_interaction):
        '''
        Compute the split potentials energies
        '''
        splitpotentials = SplitPotentialsPPI()
        try:
            splitpotentials_energies = splitpotentials.calculate_global_energies(decoy_interaction, Zscores=True)
        except RuntimeError as e:
            raise RuntimeError('Cannot compute splitpotentials energies: ' % str(e))
        return splitpotentials_energies

    def __print_energies(self, decoy_scores):
        '''
        Prints the splitpotentials and iLoops energies of each decoy
        '''
        result_fo = open(os.path.join(self.job_path, 'rmsd_Zsplitpotentials.txt'), 'w')
        for decoy in decoy_scores:
            result_fo.write('%d\t%.3f\t' % (decoy[0], decoy[1]))
            result_fo.write('%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' % (decoy[2][1]['D-PAIR'], decoy[2][1]['D-S3DC'], decoy[2][1]['D-LOCAL'], decoy[2][1]['D-3DC'], decoy[2][0]['D-3D']))
        result_fo.close()

#--------------------------#
#          MAIN            #
#--------------------------#

if __name__ == "__main__":

    options = parse_options()
   
    setup_iLoops = SetupScores(options.ppi_name, options.job_path, options.pdb_path) 
    setup_iLoops.run()