import os, sys
import numpy as np
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt

from sklearn import preprocessing
from sklearn import linear_model
from sklearn import cross_validation
from sklearn.metrics import mean_squared_error, mean_absolute_error
import cPickle

class CCharPPIAffinityParse(object):
    '''
    Get the correlation between dG and the different potentials constained in CCharPPI
    '''
    def __init__(self):
        '''
        Contructor
        '''
        self.pdb_path = os.path.expanduser('db')
        self.ba_file = os.path.join(self.pdb_path, 'affinity_benchmark1_dG.tsv')
        self.results_dir = 'R3_results'
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
        ppi_fo = open(os.path.join(self.results_dir,'ppis_used.tsv'), 'w')
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
            ppi_fo.write('{}\n'.format(ppi))
        ppi_fo.close()
        ba_fo.close()

        total_predictions_CCharPPI = {}
        total_dGs_CCharPPI = R_dGs + F_dGs
        R_predictions_CCharPPI = {}
        F_predictions_CCharPPI = {}

        enes_index = {'ZRANK':0, 'ZRANK2':1, 'ROSSETADOCK':2, 'PYDOCK':3, 'PISA':4, 'PIE':5, 'SIPPER':6}
        for x in enes_index:
            R_native_scores = [[i[enes_index[x]]] for i in R_dock_scores]
            F_native_scores = [[i[enes_index[x]]] for i in F_dock_scores]
            print '--- %s ---' % x
            print("RIGID: {}".format(len(R_dGs)))
            print("FLEXI: {}".format(len(F_dGs)))
            corr_total  = self.__traintest_linearmodel(R_native_scores, F_native_scores, R_dGs, F_dGs)

            # PREDICTION (Quim)
            # Predictions using model made from all the values
            total_predictions_CCharPPI.setdefault(x, {})
            native_scores = R_native_scores + F_native_scores
            predictions_native = self.__calculate_predictions(native_scores, corr_total)
            total_predictions_CCharPPI[x]['docking'] = '-'
            total_predictions_CCharPPI[x]['native'] = predictions_native
            # Descomposition of the predictions in Rigid and Flexible 
            R_predictions_CCharPPI.setdefault(x, {})
            F_predictions_CCharPPI.setdefault(x, {})
            for type_struc in total_predictions_CCharPPI[x]:
                R_predictions_CCharPPI[x].setdefault(type_struc, [])
                F_predictions_CCharPPI[x].setdefault(type_struc, [])
                if type_struc == 'docking':
                    R_predictions_CCharPPI[x][type_struc] = '-'
                    F_predictions_CCharPPI[x][type_struc] = '-'
                    continue
                m = 0
                for n in xrange(len(total_predictions_CCharPPI[x][type_struc])):
                    if n < len(R_dGs): # The part of the predictions list that corresponds to the length of R_dGs is RIGID
                        if total_dGs_CCharPPI[n] == R_dGs[n]:
                            R_predictions_CCharPPI[x][type_struc].append(total_predictions_CCharPPI[x][type_struc][n])
                    else: # The rest is FLEXIBLE
                        if total_dGs_CCharPPI[n] == F_dGs[m]:
                            F_predictions_CCharPPI[x][type_struc].append(total_predictions_CCharPPI[x][type_struc][n])
                            m+=1

        # Dump the predictions in a Pickle file (Quim)
        ## For the total of predictions
        dump_pred = os.path.join(self.results_dir,"predictions_CCharPPI.pcl")
        dump_dGexp = os.path.join(self.results_dir,"dGexp_CCharPPI.pcl")
        cPickle.dump(total_predictions_CCharPPI, open(dump_pred, 'w'))
        cPickle.dump(total_dGs_CCharPPI, open(dump_dGexp, 'w'))
        ## For the rigid predictions
        dump_R_pred = os.path.join(self.results_dir,"R_predictions_CCharPPI.pcl")
        dump_R_dGexp = os.path.join(self.results_dir,"R_dGexp_CCharPPI.pcl")
        cPickle.dump(R_predictions_CCharPPI, open(dump_R_pred, 'w'))
        cPickle.dump(R_dGs, open(dump_R_dGexp, 'w'))
        ## For the flexible predictions
        dump_F_pred = os.path.join(self.results_dir,"F_predictions_CCharPPI.pcl")
        dump_F_dGexp = os.path.join(self.results_dir,"F_dGexp_CCharPPI.pcl")
        cPickle.dump(F_predictions_CCharPPI, open(dump_F_pred, 'w'))
        cPickle.dump(F_dGs, open(dump_F_dGexp, 'w'))

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
        ccharppi_file = os.path.abspath('db/CCharPPI_AffinityBenchmark.csv')
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
        corr_total = [clf.coef_[0], clf.intercept_]
        return corr_total

    def __calculate_predictions(self, scores, model):
        '''
        Calculate predictions
        '''
        coef = float(model[0])
        inter = float(model[1])
        predictions = []
        for score in scores:
            if len(score) != 1: # The score is in a list, so we make sure that it is only one element before calculate the prediction
                sys.exit(10)
            pred = (coef*score[0])+inter
            predictions.append(pred)
        return predictions


#--------------------------#
#          MAIN            #
#--------------------------#

if __name__ == "__main__":
   
    ccharppi_affinityparse = CCharPPIAffinityParse() 
    ccharppi_affinityparse.run()
