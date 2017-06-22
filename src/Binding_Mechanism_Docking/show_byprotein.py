from BioLib import *
import optparse, os, sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def parse_options():
    '''
    Parse arguments
    @return options object
    '''
    parser = optparse.OptionParser('show_byprotein.py -i PPI_FILE -j JOB_PATH [-v]')
    parser.add_option("-i", dest="ppi_file", action="store", help="The ppi file [MANDATORY]", metavar="PPI_FILE")
    parser.add_option("-j", dest="job_path", action="store", help="Job path [MANDATORY]", metavar="JOB_PATH")
    parser.add_option("-v", dest="verbose", action="store_true", help="Verbose", metavar="VERBOSE")

    (options, args) = parser.parse_args()
    if options.ppi_file == None or options.job_path == None:
        parser.error('Missing arguments\n')
    return options

class ShowPPI(object):
    '''
    Print the score distributions of each location for each PPI
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
        for ppi in ppi_fo:
            ppi_type = ppi.strip('\n').split('\t')[1]
            ppi = ppi.strip('\n').split('\t')[0][0:4]
            try:
                ppi_results = self.__on_interaction(ppi)
                if self.verbose:
                    sys.stdout.write('[\033[92mSUCCESS\033[0m]\n')
                    sys.stdout.flush()
            except RuntimeError as e:
                if self.verbose:
                    sys.stdout.write('[\033[91mFAIL\033[0m] %s\n' % str(e))
                    sys.stdout.flush()
                    continue
            # Show plot
            #sns.distplot([x[4] for x in ppi_results[1]['N']], hist=False, color="g", kde_kws={"shade": True})
            sns.distplot([x[4] for x in ppi_results[1]['I']], hist=False, color="g", kde_kws={"shade": True})
            sns.distplot([x[4] for x in ppi_results[1]['P']], hist=False, color="b", kde_kws={"shade": True})
            sns.distplot([x[4] for x in ppi_results[1]['E']], hist=False, color="r", kde_kws={"shade": True})
            #plt.hist([x[0] for x in ppi_results[0]['I']], bins=100, color='g', alpha=0.5, normed=1)
            #plt.hist([x[0] for x in ppi_results[0]['P']], bins=100, color='b', alpha=0.5, normed=1)
            #plt.hist([x[0] for x in ppi_results[0]['E']], bins=100, color='r', alpha=0.5, normed=1)
            #plt.xlim(-10, 10)
            plt.show()
        ppi_fo.close()

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
            Zsp_energies = self.__get_energies(ppi, location_decoys, 'rmsd_Zsplitpotentials')
            iLoops_energies = self.__get_energies(ppi, location_decoys, 'sp_iLoops_energies')
        except Exception as e:
            raise e
        return fiberdock_energies, sp_energies, Zsp_energies, iLoops_energies

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

    def __get_size_3D(self, ppi):
        '''
        Get the size of the protein using 3D potential
        '''
        # Check if needed files exists #
        e3D_path = os.path.join(self.job_path, ppi, 'rmsd_splitpotentials.txt')
        if not os.path.isfile(e3D_path):
            raise RuntimeError('NO SPLITPOTENTIALS FILE FOUND')
        # Get 3D split potential #
        e3D_decoys = {}
        e3D_decoys_fo = open(e3D_path, 'r')
        for decoy in e3D_decoys_fo:
            decoy = decoy.strip().split('\t')
            e3D_decoys[decoy[0]] = -float(decoy[6])/1000
        e3D_decoys_fo.close()
        return e3D_decoys

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
            if etype == 'sp_iLoops_energies':
                min_e = 7
            max_e = None
            energies_path = os.path.join(self.job_path, ppi, '%s.txt' % etype)
        if not os.path.isfile(energies_path):
            raise RuntimeError('NO %s ENERGIES FILE FOUND' % etype.upper())
        e3D = self.__get_size_3D(ppi)
        # Get energies #
        energies = {'N': [], 'I': [], 'P': [], 'E': []}
        if etype == 'sp_iLoops_energies':
            energies = {'N': [[np.nan]*6], 'I': [[np.nan]*6], 'P': [[np.nan]*6], 'E': [[np.nan]*6]}
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
        return energies

#--------------------------#
#          MAIN            #
#--------------------------#

if __name__ == "__main__":

    options = parse_options()

    ppi_file = os.path.abspath(options.ppi_file)
    job_path = os.path.abspath(options.job_path)
   
    show_ppi = ShowPPI(ppi_file, job_path, options.verbose) 
    show_ppi.run()
