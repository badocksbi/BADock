from BioLib import *
import optparse, os, sys, random
import numpy as np

def parse_options():
    '''
    Parse arguments
    @return options object
    '''
    parser = optparse.OptionParser('analyze_byprotein.py -i PPI_FILE -j JOB_PATH [-v]')
    parser.add_option("-i", dest="ppi_file", action="store", help="The ppi file [MANDATORY]", metavar="PPI_FILE")
    parser.add_option("-j", dest="job_path", action="store", help="Job path [MANDATORY]", metavar="JOB_PATH")
    parser.add_option("-v", dest="verbose", action="store_true", help="Verbose", metavar="VERBOSE")

    (options, args) = parser.parse_args()
    if options.ppi_file == None or options.job_path == None:
        parser.error('Missing arguments\n')
    return options

class AnalyzePPI(object):
    '''
    Get the scores of each decoy for each protein-protein interaction and analyze the
    diferent locations
    '''
    def __init__(self, ppi_file, job_path, verbose):
        '''
        Contructor
        '''
        self.ppi_file = ppi_file
        self.job_path = job_path
        self.verbose = verbose
        
    def run(self):
        '''
        Runs the computations
        '''
        ppi_fo = open(self.ppi_file, 'r')
        energies_print = []
        for ppi in ppi_fo:
            ppi_type = ppi.strip('\n').split('\t')[1]
            if ppi_type[0] != 'E' and ppi_type != 'OX':
                continue
            ppi = ppi.strip('\n').split('\t')[0][0:4]
            try:
                ppi_results = self.__on_interaction(ppi)
                energies_print.append(ppi_results)
                if self.verbose:
                    sys.stdout.write('[\033[92mSUCCESS\033[0m]\n')
                    sys.stdout.flush()
            except RuntimeError as e:
                if self.verbose:
                    sys.stdout.write('[\033[91mFAIL\033[0m] %s\n' % str(e))
                    sys.stdout.flush()
        ppi_fo.close()
        self.__print(energies_print)

    def __on_interaction(self, ppi):
        '''
        Computations for each ppi
        '''
        if self.verbose:
            sys.stdout.write('PPI %s --> ' % ppi)
            sys.stdout.flush()
        # Get the energies by location #
        try:
            location_decoys = self.__get_location_decoys(ppi)
        except RuntimeError as e:
            raise e
        if not location_decoys['N'] or not location_decoys['I'] or not location_decoys['P'] or not location_decoys['E']:
            raise RuntimeError('NO NATIVE, INTERFACE, PARTIAL OR EXTERNAL DECOYS')
        try:
            fiberdock_energies = self.__get_energies(ppi, location_decoys, 'fiberdock')
            sp_energies = self.__get_energies(ppi, location_decoys, 'rmsd_splitpotentials')
        except Exception as e:
            raise e
        # Compute average macrostate energies #
        macrostate_N = np.concatenate((np.mean(fiberdock_energies['N'], axis=0), 
                                       np.mean(sp_energies['N'], axis=0),
                                       np.std(fiberdock_energies['N'], axis=0),
                                       np.std(sp_energies['N'], axis=0)))
        macrostate_I = np.concatenate((np.mean(fiberdock_energies['I'], axis=0), 
                                       np.mean(sp_energies['I'], axis=0),
                                       np.std(fiberdock_energies['I'], axis=0),
                                       np.std(sp_energies['I'], axis=0)))
        macrostate_P = np.concatenate((np.mean(fiberdock_energies['P'], axis=0), 
                                       np.mean(sp_energies['P'], axis=0),
                                       np.std(fiberdock_energies['P'], axis=0),
                                       np.std(sp_energies['P'], axis=0)))
        macrostate_E = np.concatenate((np.mean(fiberdock_energies['E'], axis=0), 
                                       np.mean(sp_energies['E'], axis=0),
                                       np.std(fiberdock_energies['E'], axis=0),
                                       np.std(sp_energies['E'], axis=0)))
        return ppi, macrostate_N, macrostate_I, macrostate_P, macrostate_E

    def __get_rmsd(self, ppi):
        '''
        Get the RMSD for ach decoy using the sp_iLoops_energies file
        '''
        # Check if needed files exists #
        rmsd_path = os.path.join(self.job_path, ppi, 'rmsd_splitpotentials.txt')
        if not os.path.isfile(rmsd_path):
            raise RuntimeError('NO RMSD FILE FOUND')
        # Get rmsd #
        rmsd_decoys = {}
        rmsd_decoys_fo = open(rmsd_path, 'r')
        for decoy in rmsd_decoys_fo:
            decoy = decoy.strip().split('\t')
            rmsd_decoys[decoy[0]] = float(decoy[1])
        rmsd_decoys_fo.close()
        return rmsd_decoys

    def __get_location_decoys(self, ppi):
        '''
        Get the decoys for each location
        '''
        # Check if needed files exists #
        decoys_location_path = os.path.join(self.job_path, ppi, 'decoys_location.txt')
        if not os.path.isfile(decoys_location_path):
            raise RuntimeError('NO DECOYS LOCATION FILE FOUND')
        # Get rmsds #
        location_decoys = {'N': [],'I': [], 'P': [], 'E': []}
        try:
            rmsd_decoys = self.__get_rmsd(ppi)
        except RuntimeError as e:
            raise e
        # Get location decoys #
        decoys_locations_fo = open(decoys_location_path, 'r')
        for decoy in decoys_locations_fo:
            decoy = decoy.strip('\n').split('\t')
            if not decoy[0] in rmsd_decoys:
                continue
            if rmsd_decoys[decoy[0]] < 10:
                location_decoys['N'].append(decoy[0])
            else:
                location_decoys[decoy[1]].append(decoy[0])
        decoys_locations_fo.close()
        return location_decoys

    def __get_energies(self, ppi, location_decoys, etype, zscores=False):
        '''
        Gets the energies for each location
        etype: fiberdock, sp or iLoops for choosing the energie type
        '''
        if etype == 'fiberdock':
            separator = '|'
            min_e = 1
            max_e = 11
            energies_path = os.path.join(self.job_path, ppi, '%s_energies.txt' % etype)
            energies_path_ref = os.path.join(self.job_path, ppi, '%s_energies.ref' % etype)
            if os.path.isfile(energies_path_ref):
                energies_path = energies_path_ref
        else:
            separator = '\t'
            min_e = 2
            max_e = None
            energies_path = os.path.join(self.job_path, ppi, '%s.txt' % etype)
        if not os.path.isfile(energies_path):
            raise RuntimeError('NO %s ENERGIES FILE FOUND' % etype.upper())
        # Get energies #
        energies = {'N': [], 'I': [], 'P': [], 'E': []}
        energies_fo = open(energies_path, 'r')
        for decoy in energies_fo:
            decoy = decoy.strip('\n').split(separator)
            if decoy[0].strip() in location_decoys['N']:
                energies['N'].append([float(energy.strip()) for energy in decoy[min_e:max_e]])
            if decoy[0].strip() in location_decoys['I']:
                energies['I'].append([float(energy.strip()) for energy in decoy[min_e:max_e]])
            if decoy[0].strip() in location_decoys['P']:
                energies['P'].append([float(energy.strip()) for energy in decoy[min_e:max_e]])
            if decoy[0].strip() in location_decoys['E']:
                energies['E'].append([float(energy.strip()) for energy in decoy[min_e:max_e]])
        energies_fo.close()
        if zscores:
            all_energies = energies['N']+energies['I']+energies['P']+energies['E']
            zenergies = {'N': [], 'I': [], 'P': [], 'E': []}
            zmean = np.nanmean(all_energies, axis=0)
            zstd = np.nanstd(all_energies, axis=0)
            for loc in energies:
                for decoy in energies[loc]:
                    zenergies[loc].append((decoy-zmean)/zstd)
            return zenergies
        else:
            return energies

    def __print(self, energies):
        '''
        Print the resultant energies
        '''
        ene_fo = open('scores_locations_E+OX_std.txt' , 'w')
        for energie in energies:
            ene_fo.write('%s\t' % energie[0])
            for x in xrange(len(energie[1])):
                ene_fo.write('%.2f\t%.2f\t%.2f\t%.2f' % (energie[1][x], energie[2][x], energie[3][x], energie[4][x]))
                if x < len(energie[1])-1:
                    ene_fo.write('\t')
            ene_fo.write('\n')
        ene_fo.close()

#--------------------------#
#          MAIN            #
#--------------------------#

if __name__ == "__main__":

    options = parse_options()

    ppi_file = os.path.abspath(options.ppi_file)
    job_path = os.path.abspath(options.job_path)
   
    analyze_ppi = AnalyzePPI(ppi_file, job_path, options.verbose) 
    analyze_ppi.run()
