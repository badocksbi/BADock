from BioLib import *
import os, sys

class ScoreNative(object):
    '''
    Score the native bound complexes
    '''
    def __init__(self):
        '''
        Contructor
        '''
        self.pdb_path = os.path.expanduser('~/Databases/AffinityBenchmark')
        self.ba_file = os.path.join(self.pdb_path, 'affinity_benchmark_all_dG.tsv')
        self.verbose = True

    def run(self):
        '''
        Runs the computations
        '''
        ba_fo = open(self.ba_file, 'r')
        scores_native_fo = open('scores_native.txt', 'w')
        for ppi in ba_fo:
            ppi = ppi.strip('\n').split('\t')
            ppi = ppi[0]
            if self.verbose:
                sys.stdout.write('PPI %s --> ' % ppi)
                sys.stdout.flush()
            try:
                dock_score = self.__on_interaction(ppi)
            except Exception as e:
                if self.verbose:
                    sys.stdout.write('FAIL: %s\n' % str(e))
                    sys.stdout.flush()
                continue
            scores_str = '%s\t' % ppi
            scores_str += '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t' % tuple(dock_score[0])
            scores_str += '%d\t%d\t%d\t%d\t%d\t%d\t' % (dock_score[1]['charged-charged'], dock_score[1]['charged-polar'], dock_score[1]['charged-apolar'], 
                                                        dock_score[1]['polar-polar'], dock_score[1]['polar-apolar'], dock_score[1]['apolar-apolar'])
            scores_str += '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (dock_score[2][0]['D-PAIR'], dock_score[2][0]['D-S3DC'], dock_score[2][0]['D-LOCAL'], dock_score[2][0]['D-3DC'], dock_score[2][0]['D-3D'])
            scores_native_fo.write(scores_str)
            scores_native_fo.flush()
            if self.verbose:
                sys.stdout.write('SUCCESS\n')
                sys.stdout.flush()
        scores_native_fo.close()
        ba_fo.close()

    def __on_interaction(self, ppi):
        '''
        Computations for each protein-protein interaction
        '''        
        # GET THE STRUCTURES FROM THE PDB FILES #
        receptor_pdb = os.path.join(self.pdb_path, ppi+'_r_b.pdb')
        ligand_pdb = os.path.join(self.pdb_path, ppi+'_l_b.pdb')
        try:
            struct_rec, struct_lig = self.__get_structures(receptor_pdb, ligand_pdb)
        except RuntimeError as e:
            raise e
        # GET THE NATIVE SCORES #
        native_interaction = Interaction(struct_rec, struct_lig)
        try:
            fiberdock_energies = self.__get_fiberdock(ppi, receptor_pdb, ligand_pdb)
            ic_scores = IC_scores(native_interaction).calculate()
            splitpotentials_energies = self.__get_splitpotentials(native_interaction)
        except RuntimeError as e:
            raise e
        return fiberdock_energies, ic_scores, splitpotentials_energies

    def __get_structures(self, receptor_pdb, ligand_pdb):
        '''
        Parse the original PDB files in order to get the structure objects
        '''
        try:
            struct_rec = PDB.read_pdb(receptor_pdb, merge_chains=True)
            struct_rec.clean()
            struct_rec.set_dssp()
            struct_lig = PDB.read_pdb(ligand_pdb, merge_chains=True)
            struct_lig.clean()
            struct_lig.set_dssp()
        except Exception as e:
            raise RuntimeError('ERROR: Cannot parse PDB: %s\n' % str(e))
        return struct_rec, struct_lig

    def __get_splitpotentials(self, interaction):
        '''
        Compute the split potentials energies
        '''
        splitpotentials = SplitPotentialsPPI()
        try:
            splitpotentials_energies = splitpotentials.calculate_global_energies(interaction, Zscores=False)
        except RuntimeError as e:
            raise RuntimeError('Cannot compute splitpotentials energies: ' % str(e))
        return splitpotentials_energies

    def __get_fiberdock(self, ppi, receptor_pdb, ligand_pdb):
        '''
        Compute the fiberdock energies
        '''
        # Execute FiberDock #
        job_path = 'ba_jobs_u/%s' % ppi
        self.__create_fiberdock_patchdock(job_path)
        fiberdock_params_path = os.path.join(job_path, 'fiberdock_params_native.txt')
        fiberdock_params_fo = open(fiberdock_params_path, 'w')
        fiberdock_params_fo.write(self.__create_fiberdock_params(job_path, receptor_pdb, ligand_pdb))
        fiberdock_params_fo.close()
        fiberdock_exe_cmd = '/soft/bio/structure/FiberDock1.1/FiberDock %s' % fiberdock_params_path
        try:
            subprocess.check_output(fiberdock_exe_cmd, stderr=subprocess.STDOUT, shell=True)
        except Exception as e:
            raise RuntimeError('ERROR: Error executing FiberDock: %s\n' % str(e))
        os.remove(fiberdock_params_path)
        fiberdock_energies_fo = open(os.path.join(job_path, 'fiberdock_energies_native.ref'))
        for native in fiberdock_energies_fo:
            native = native.strip('\n').split('|')
            if native[0].strip() == '1':
                return [float(energy.strip()) for energy in native[1:11]]

    def __create_fiberdock_patchdock(self, job_path):
        '''
        '''
        patchdock_native_path = os.path.join(job_path, 'patchdock_native.txt')
        patchdock_native_fo = open(patchdock_native_path, 'w')
        patchdock_native_fo.write('1 0 0 0 0 0 0\n')
        patchdock_native_fo.close()

    def __create_fiberdock_params(self, job_path, receptor_pdb, ligand_pdb):
        '''
        Returns the fiberdock params string
        '''
        ### I/O
        params =  'receptorPDBFileName %s.HB\n' %  receptor_pdb
        params += 'ligandPDBFileName %s.HB\n' % ligand_pdb
        # transformations for refinement
        params += 'transFileName %s\n' % os.path.join(job_path, 'patchdock_native.txt')
        ## libraries files
        params += 'rotamerLibFile  /soft/bio/structure/FiberDock1.1/lib/bbdep02.May.sortlib\n'
        params += 'protLib /soft/bio/structure/FiberDock1.1/lib/chem.lib\n'
        ## output file
        params += 'energiesOutFileName %s\n' % os.path.join(job_path, 'fiberdock_energies_native')
        ### Output Options
        # to output refined complexes
        params += 'printRefinedComplexes 0\n'
        ## 1 - only energy caclulaltion is performed without refinement (works only for FiberDock pre)
        params += 'onlyEnergyCalculation 1\n'
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
        params += 'rigidBodyMinimizationCycles 0\n'
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
        params += 'performBackboneRef 0\n'
        #params += 'Cycles 10\n'
        #params += 'RecAnalyzedFile %s.CA\n' % receptor_pdb
        #params += 'ReceptorNMs %s_r_b.nma\n' % os.path.join(job_path, ppi)
        #params += 'Lambda 0.05\n'
        #params += 'LigAnalyzedFile %s.CA\n' % ligand_pdb
        #params += 'LigandNMs %s_l_b.nma\n' % os.path.join(job_path, ppi)
        ## LOG-FILE params
        params += 'log-file %s\n' % os.path.join(job_path, 'fiberdock_log_native.txt')
        params += 'log-level 2\n'
        return params

#--------------------------#
#          MAIN            #
#--------------------------#

if __name__ == "__main__":
   
    Score_native = ScoreNative() 
    Score_native.run()