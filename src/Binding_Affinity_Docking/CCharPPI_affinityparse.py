import os, sys
import numpy as np
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt

from sklearn import preprocessing
from sklearn import linear_model
from sklearn import cross_validation
from sklearn.metrics import mean_squared_error, mean_absolute_error

class CCharPPIAffinityParse(object):
    '''
    Get the correlation between dG and the different potentials constained in CCharPPI
    '''
    def __init__(self):
        '''
        Contructor
        '''
        self.pdb_path = os.path.expanduser('~/Databases/AffinityBenchmark')
        self.ba_file = os.path.join(self.pdb_path, 'affinity_benchmark1_dG.tsv')
        self.verbose = True

    def run(self):
        '''
        Runs the computations
        '''
        ccharppi_scores = self.__get_CCharPPI()
        R_dock_scores = []
        R_dGs = []
        F_dock_scores = []
        F_dGs = []
        ba_fo = open(self.ba_file, 'r')
        for ppi in ba_fo:
            ppi = ppi.strip('\n').split('\t')
            ppi_type = ppi[1]
            dG = float(ppi[3])
            irmsd = float(ppi[4])
            ppi = ppi[0]
            if ppi_type[0] != 'E' and ppi_type != 'OX':
                continue
            #if irmsd > 1.0:
            #    continue
            if self.verbose:
                sys.stdout.write('PPI %s --> ' % ppi)
                sys.stdout.flush()
            if ppi not in ccharppi_scores:
                if self.verbose:
                    sys.stdout.write('FAIL: NOT FOUND IN CCHARPPI FILE\n')
                    sys.stdout.flush()
                continue
            if irmsd <= 1.0:
                R_dock_scores.append(ccharppi_scores[ppi])
                R_dGs.append(dG)
            else:
                F_dock_scores.append(ccharppi_scores[ppi])
                F_dGs.append(dG)
            if self.verbose:
                sys.stdout.write('SUCCESS\n')
                sys.stdout.flush()
        ba_fo.close()

        enes_index = {'ZRANK':0, 'ZRANK2':1, 'ROSSETADOCK':2, 'PYDOCK':3, 'PISA':4, 'PIE':5, 'SIPPER':6}
        for x in enes_index:
            R_native_scores = [[i[enes_index[x]]] for i in R_dock_scores]
            F_native_scores = [[i[enes_index[x]]] for i in F_dock_scores]
            print '--- %s ---' % x
            self.__traintest_linearmodel(R_native_scores, F_native_scores, R_dGs, F_dGs)
        exit()


        sys.stdout.write('# Complexes: %d\n' % len(dGs))
        sys.stdout.write('ZRANK %.2f\t%.2e\n' % pearsonr([i[0] for i in dock_scores], dGs))
        sys.stdout.write('ZRANK2 %.2f\t%.2e\n' % pearsonr([i[1] for i in dock_scores], dGs))
        sys.stdout.write('ROSSETADOCK %.2f\t%.2e\n' % pearsonr([i[2] for i in dock_scores], dGs))
        sys.stdout.write('PYDOCK %.2f\t%.2e\n' % pearsonr([i[3] for i in dock_scores], dGs))
        sys.stdout.write('PISA %.2f\t%.2e\n' % pearsonr([i[4] for i in dock_scores], dGs))
        sys.stdout.write('PIE %.2f\t%.2e\n' % pearsonr([i[5] for i in dock_scores], dGs))
        sys.stdout.write('SIPPER %.2f\t%.2e\n' % pearsonr([i[6] for i in dock_scores], dGs))
        #sys.stdout.write('FIREDOCK %.2f\t%.2e\n' % pearsonr([i[7] for i in dock_scores], dGs))

    def __get_CCharPPI(self):
        '''
        Get the CCharPPI potentials
        '''
        ppis = {}
        ccharppi_file = os.path.abspath('CCharPPI_AffinityBenchmark.csv')
        ccharppi_fo = open(ccharppi_file, 'r')
        ccharppi_fo.readline()
        for ppi in ccharppi_fo:
            ppi = ppi.strip('\n').split(',')
            ppis[ppi[0][0:4]] = (float(ppi[47]), float(ppi[48]), float(ppi[62]), float(ppi[67]), float(ppi[101]), float(ppi[105]), float(ppi[70]), float(ppi[102]))
        ccharppi_fo.close()
        return ppis

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
            Fpearsons.append(pearsonr(Fy_predict, Fy_test)[0])
            Fmses.append(mean_squared_error(Fy_predict, Fy_test)**0.5)
        print 'TOTAL -> Pearsonr: %.2f (std %.2f), RMSE: %.2f' % (np.nanmean(Tpearsons), np.nanstd(Tpearsons), np.nanmean(Tmses))
        print 'RIGID -> Pearsonr: %.2f (std %.2f), RMSE: %.2f' % (np.nanmean(Rpearsons), np.nanstd(Rpearsons), np.nanmean(Rmses))
        print 'FLEXI -> Pearsonr: %.2f (std %.2f), RMSE: %.2f' % (np.nanmean(Fpearsons), np.nanstd(Fpearsons), np.nanmean(Fmses))
        # Create model
        clf = linear_model.LinearRegression()
        clf.fit(R_scores+F_scores, R_dG+F_dG)
        print 'Coefs: %.10f*E, Intercept: %.3f' % (clf.coef_[0], clf.intercept_)

#--------------------------#
#          MAIN            #
#--------------------------#

if __name__ == "__main__":
   
    ccharppi_affinityparse = CCharPPIAffinityParse() 
    ccharppi_affinityparse.run()
