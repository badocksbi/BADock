import os, sys
import numpy as np
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import preprocessing
from sklearn import linear_model
from sklearn import cross_validation
from sklearn.metrics import mean_squared_error, mean_absolute_error

class AffinityAnalysis(object):
    '''
    Analyze the correlation between the binding affinities and the different parameters
    computed for each ppi docking or complex
    '''
    def __init__(self):
        '''
        Contructor
        '''
        self.jobs_path = os.path.abspath('ba_jobs')
        self.pdb_path = os.path.expanduser('~/Databases/AffinityBenchmark')
        self.ba_file = os.path.join(self.pdb_path, 'affinity_benchmark_all_dG.tsv')
        self.analysis = 'Native' # 'Native' or 'Docking'
        self.filter_type = 'D' # 'E', 'A', 'O' or None for all or 'D' for E+OX (dataset)
        self.filter_irmds = None # 'R', 'F' or None for all
        self.calculations = None # 'KON' 'KOFF' or None for dGs or Kds

    def run(self):
        '''
        Runs the computations
        '''
        native_scores = []
        docking_scores = []
        dGs = []
        kds = []
        kons = []
        koffs = []
        R_native = []
        R_docking = []
        F_native = []
        F_docking = []
        R_dG = []
        F_dG = []
        ba_fo = open(self.ba_file, 'r')
        for ppi in ba_fo:
            # Get ppf features # 
            ppi = ppi.strip('\n').split('\t')
            ppi_type = ppi[1]
            irmsd = float(ppi[4])
            kd = np.log(float(ppi[2]))
            dG = float(ppi[3])
            kon = -np.log(float(ppi[5]))
            koff = np.log(float(ppi[6]))
            # Filter by ppi type and irmsd and nan values for kon-koff #
            if self.calculations == 'KON' and np.isnan(kon):
                continue
            if self.calculations == 'KOFF' and np.isnan(koff):
                continue
            ppi = ppi[0]
            if self.filter_type and self.filter_type != 'D' and ppi_type[0] != self.filter_type:
                continue
            if self.filter_type == 'D' and ppi_type[0] != 'E' and ppi_type != 'OX':
                continue
            if self.filter_irmds == 'R' and irmsd > 1.0:
                continue
            if self.filter_irmds == 'F' and irmsd < 1.0:
                continue
            # End filtering #
            #sys.stdout.write('PPI %s --> ' % ppi)
            #sys.stdout.flush()
            scores_native_fo = open('scores_native.txt')
            for interaction in scores_native_fo:
                interaction = interaction.strip('\n').split('\t')
                if interaction[0] == ppi:
                    scores_native = [float(energy) for energy in interaction[1:]]
                    break
            scores_native_fo.close()
            try:
                scores_docking = self.__on_docking(ppi)
            except Exception as e:
                sys.stdout.write('FAIL: %s\n' % str(e))
                sys.stdout.flush()
                continue
            native_scores.append(scores_native)
            docking_scores.append(scores_docking)
            if irmsd <= 1.0:
                R_native.append(scores_native)
                R_docking.append(scores_docking)
                R_dG.append(dG)
            if irmsd > 1.0:
                F_native.append(scores_native)
                F_docking.append(scores_docking)
                F_dG.append(dG)
            dGs.append(dG)
            kds.append(kd)
            kons.append(kon)
            koffs.append(koff)
            #sys.stdout.write('SUCCESS\n')
            #sys.stdout.flush()
        ba_fo.close()

        #print pearsonr([i[5] for i in native_scores], [i[5] for i in docking_scores])
        if self.analysis == 'Native':
            ppi_scores = native_scores
        else:
            ppi_scores = docking_scores
        if self.calculations == 'KON':
            self.__print_results(ppi_scores, kons, R_native, R_docking, F_native, F_docking, R_dG, F_dG)
        elif self.calculations == 'KOFF':
            self.__print_results(ppi_scores, koffs, R_native, R_docking, F_native, F_docking, R_dG, F_dG)
        else:
            self.__print_results(ppi_scores, dGs, R_native, R_docking, F_native, F_docking, R_dG, F_dG)

    def __on_docking(self, ppi):
        '''
        Computations for each protein-protein interaction
        '''
        scores_fo = open(os.path.join(self.jobs_path, ppi, 'decoys_scores.txt'), 'r')
        scores = []
        for decoy in scores_fo:
            decoy = [float(x) for x in decoy.strip().split('\t')[2:]]
            scores.append(decoy)
        scores_fo.close()
        scores_mean = list(np.mean(scores, axis=0))
        #scores_std = list(np.std(scores, axis=0))
        return scores_mean

    def __print_results(self, dock_scores, dGs, R_native, R_docking, F_native, F_docking, R_dG, F_dG):
        '''
        Print correlations and predictive models
        '''
        # Fiberdock #
        sys.stdout.write('FiberDock %.2f\t%.2e\n' % pearsonr([i[0] for i in dock_scores], dGs))
        sys.stdout.write('aVdW %.2f\t%.2e\n' % pearsonr([i[1] for i in dock_scores], dGs))
        sys.stdout.write('rVdW %.2f\t%.2e\n' % pearsonr([i[2] for i in dock_scores], dGs))
        #sys.stdout.write('ACE %.2f\t%.2e\n' % pearsonr([i[3] for i in dock_scores], dGs))
        #sys.stdout.write('inside %.2f\t%.2e\n' % pearsonr([i[4] for i in dock_scores], dGs))
        sys.stdout.write('aElec %.2f\t%.2e\n' % pearsonr([i[5] for i in dock_scores], dGs))
        sys.stdout.write('rElec %.2f\t%.2e\n' % pearsonr([i[6] for i in dock_scores], dGs))
        sys.stdout.write('laElec %.2f\t%.2e\n' % pearsonr([i[7] for i in dock_scores], dGs))
        sys.stdout.write('lrElec %.2f\t%.2e\n' % pearsonr([i[8] for i in dock_scores], dGs))
        sys.stdout.write('HB %.2f\t%.2e\n' % pearsonr([i[9] for i in dock_scores], dGs))
        # ICs #
        #sys.stdout.write('Total %.2f\t%.2e\n' % pearsonr([i[10]+i[11]+i[12]+i[13]+i[14]+i[15] for i in dock_scores], dGs))
        #sys.stdout.write('Charged-Charged %.2f\t%.2e\n' % pearsonr([i[10] for i in dock_scores], dGs))
        #sys.stdout.write('Charged-Polar %.2f\t%.2e\n' % pearsonr([i[11] for i in dock_scores], dGs))
        #sys.stdout.write('Charged-Apolar %.2f\t%.2e\n' % pearsonr([i[12] for i in dock_scores], dGs))
        #sys.stdout.write('Polar-Polar %.2f\t%.2e\n' % pearsonr([i[13] for i in dock_scores], dGs))
        #sys.stdout.write('Polar-Apolar %.2f\t%.2e\n' % pearsonr([i[14] for i in dock_scores], dGs))
        #sys.stdout.write('Apolar-Apolar %.2f\t%.2e\n' % pearsonr([i[15] for i in dock_scores], dGs))
        # Split potentials #
        sys.stdout.write('EPAIR %.2f\t%.2e\n' % pearsonr([i[16] for i in dock_scores], dGs))
        sys.stdout.write('ES3DC %.2f\t%.2e\n' % pearsonr([i[17] for i in dock_scores], dGs))
        #sys.stdout.write('ELOCAL %.2f\t%.2e\n' % pearsonr([i[18] for i in dock_scores], dGs))
        #sys.stdout.write('E3DC %.2f\t%.2e\n' % pearsonr([i[19] for i in dock_scores], dGs))
        sys.stdout.write('E3D %.2f\t%.2e\n' % pearsonr([i[20] for i in dock_scores], dGs))
        #print len(dGs)

        enes_index = {'FiberDock':0, 'aVdW':1, 'rVdW':2, 'aElec':5, 'rElec':6, 'laElec':7, 'lrElec':8, 'HB':9, 'EPAIR':16, 'ES3DC':17, 'E3D':20}
        for x in enes_index:
            R_native_scores = [[i[enes_index[x]]] for i in R_native]
            F_native_scores = [[i[enes_index[x]]] for i in F_native]
            R_docking_scores = [[i[enes_index[x]]] for i in R_docking]
            F_docking_scores = [[i[enes_index[x]]] for i in F_docking]
            print '--- %s --- NATIVE ---' % x
            self.__traintest_linearmodel(R_native_scores, F_native_scores, R_dG, F_dG)
            print '--- %s --- DOCKING --' % x
            self.__traintest_linearmodel(R_docking_scores, F_docking_scores, R_dG, F_dG)

        R_docking_scores = [[i[enes_index['ES3DC']]] for i in R_docking]
        F_docking_scores = [[i[enes_index['ES3DC']]] for i in F_docking]
        R_predicted_dG, R_experimental_dG, F_predicted_dG, F_experimental_dG = self.__traintest_linearmodel(R_docking_scores, F_docking_scores, R_dG, F_dG)

        sns.set_context("paper", font_scale=1.5)
        sns.set_style("white")
        plt.ylim(-4,-19)
        plt.xlim(-8, -13.5)
        plt.xlabel(r'Predicted $\Delta$Gs (kcal mol$^{-1}$)')
        plt.ylabel(r'Experimental $\Delta$Gs (kcal mol$^{-1}$)')

        #sns.jointplot(np.array(R_predicted_dG+F_predicted_dG), np.array(R_experimental_dG+F_experimental_dG), kind="kde");
        sns.kdeplot(np.array(R_predicted_dG), np.array(R_experimental_dG), shade_lowest=False, cmap="Blues")
        sns.kdeplot(np.array(F_predicted_dG), np.array(F_experimental_dG), shade_lowest=False, cmap="Reds")
        #plt.scatter(R_predicted_dG, R_experimental_dG, color='b', marker='o', label='Rigid')
        #plt.scatter(F_predicted_dG, F_experimental_dG, color='r', marker='v', label='Flexible')
        #legend = plt.legend(loc='lower right', frameon=True)
        #legend.get_frame().set_facecolor('#E1E1E1')
        plt.savefig('img/BA_ES3DC_Docking.eps', bbox_inches='tight', dpi=350)
        plt.show()
        
    def __traintest_linearmodel(self, R_scores, F_scores, R_dG, F_dG):
        '''
        Print linear model stats
        '''
        Tpearsons = []
        Tmses = []
        Rpearsons = []
        Rmses = []
        Fpearsons = []
        Fmses = []
        R_predicted_dG = []
        R_experimental_dG = []
        F_predicted_dG = []
        F_experimental_dG = []       
        for x in xrange(1000):
            clf = linear_model.LinearRegression()
            RX_train, RX_test, Ry_train, Ry_test = cross_validation.train_test_split(R_scores, R_dG, test_size=0.10, random_state=x+42)
            FX_train, FX_test, Fy_train, Fy_test = cross_validation.train_test_split(F_scores, F_dG, test_size=0.10, random_state=x+42)
            clf.fit(RX_train+FX_train, Ry_train+Fy_train)
            Ty_predict = clf.predict(RX_test+FX_test)
            Ry_predict = clf.predict(RX_test)
            Fy_predict = clf.predict(FX_test)
            Tpearsons.append(pearsonr(Ty_predict, Ry_test+Fy_test)[0])
            Tmses.append(mean_squared_error(Ty_predict, Ry_test+Fy_test)**0.5)
            Rpearsons.append(pearsonr(Ry_predict, Ry_test)[0])
            Rmses.append(mean_squared_error(Ry_predict, Ry_test)**0.5)
            R_predicted_dG.extend(Ry_predict)
            R_experimental_dG.extend(Ry_test)
            Fpearsons.append(pearsonr(Fy_predict, Fy_test)[0])
            Fmses.append(mean_squared_error(Fy_predict, Fy_test)**0.5)
            F_predicted_dG.extend(Fy_predict)
            F_experimental_dG.extend(Fy_test)     
        print 'TOTAL -> Pearsonr: %.2f (std %.2f), RMSE: %.2f' % (np.nanmean(Tpearsons), np.nanstd(Tpearsons), np.nanmean(Tmses))
        print 'RIGID -> Pearsonr: %.2f (std %.2f), RMSE: %.2f' % (np.nanmean(Rpearsons), np.nanstd(Rpearsons), np.nanmean(Rmses))
        print 'FLEXI -> Pearsonr: %.2f (std %.2f), RMSE: %.2f' % (np.nanmean(Fpearsons), np.nanstd(Fpearsons), np.nanmean(Fmses))
        # Create model
        clf = linear_model.LinearRegression()
        clf.fit(R_scores+F_scores, R_dG+F_dG)
        print 'Coefs: %.10f*E, Intercept: %.3f' % (clf.coef_[0], clf.intercept_)
        return R_predicted_dG, R_experimental_dG, F_predicted_dG, F_experimental_dG

#--------------------------#
#          MAIN            #
#--------------------------#

if __name__ == "__main__":
   
    affinity_analysis = AffinityAnalysis() 
    affinity_analysis.run()
