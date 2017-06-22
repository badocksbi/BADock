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
        self.jobs_path = os.path.abspath('/Users/quim/Documents/ba_jobs')
        self.pdb_path = os.path.expanduser('db')
        self.ba_file = os.path.join(self.pdb_path, 'affinity_benchmark_all_dG.tsv')
        self.analysis = 'Native' # 'Native' or 'Docking'
        self.filter_type = 'D' # 'E', 'A', 'O' or None for all or 'D' for E+OX (dataset)
        self.filter_irmds = None # 'R', 'F' or None for all
        self.calculations = None # 'KON' 'KOFF' or None for dGs or Kds
        self.img_path = 'img' # Path were you want to save the resulting pictures (end the path without '/')
        self.distance_reg_line = 1.4 # This is the distance to the regression line of the lines that will define the area of best regression


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
            scores_native_fo = open('db/scores_native.txt')
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

        percents_doc_file = 'percents_all_decoys.tsv'
        percents_nat_file = 'percents_native.tsv'
        fpd = open(percents_doc_file, 'w')
        fpn = open(percents_nat_file, 'w')

        absolute_doc_file = 'absolute_all_decoys.tsv'
        absolute_nat_file = 'absolute_native.tsv'
        fad = open(absolute_doc_file, 'w')
        fan = open(absolute_nat_file, 'w')

        enes_index = {'FiberDock':0, 'aVdW':1, 'rVdW':2, 'aElec':5, 'rElec':6, 'laElec':7, 'lrElec':8, 'HB':9, 'EPAIR':16, 'ES3DC':17, 'E3D':20}
        # for x in enes_index:
        #     R_native_scores = [[i[enes_index[x]]] for i in R_native]
        #     F_native_scores = [[i[enes_index[x]]] for i in F_native]
        #     R_docking_scores = [[i[enes_index[x]]] for i in R_docking]
        #     F_docking_scores = [[i[enes_index[x]]] for i in F_docking]
        #     print '--- %s --- NATIVE ---' % x
        #     R_predicted_dG_nat, R_experimental_dG_nat, F_predicted_dG_nat, F_experimental_dG_nat, T_filt_nat, R_filt_nat, F_filt_nat, percents_nat, absolute_nat, _ = self.__traintest_linearmodel(R_native_scores, F_native_scores, R_dG, F_dG)
        #     print '--- %s --- DOCKING --' % x
        #     R_predicted_dG, R_experimental_dG, F_predicted_dG, F_experimental_dG, T_filt, R_filt, F_filt, percents_doc, absolute_doc, _ = self.__traintest_linearmodel(R_docking_scores, F_docking_scores, R_dG, F_dG)

        #     # Quim: Writing an output file with the percentages of differences dGexp - dGpred (Native)
        #     fpn.write('%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' %( x,
        #                         percents_nat['T']['1_4'], percents_nat['T']['2_8'], percents_nat['T']['4_2'], percents_nat['T']['more'],
        #                         percents_nat['R']['1_4'], percents_nat['R']['2_8'], percents_nat['R']['4_2'], percents_nat['R']['more'],
        #                         percents_nat['F']['1_4'], percents_nat['F']['2_8'], percents_nat['F']['4_2'], percents_nat['F']['more'] ))

        #     # Quim: Writing an output file with the percentages of differences dGexp - dGpred (All decoys)
        #     fpd.write('%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' %( x,
        #                         percents_doc['T']['1_4'], percents_doc['T']['2_8'], percents_doc['T']['4_2'], percents_doc['T']['more'],
        #                         percents_doc['R']['1_4'], percents_doc['R']['2_8'], percents_doc['R']['4_2'], percents_doc['R']['more'],
        #                         percents_doc['F']['1_4'], percents_doc['F']['2_8'], percents_doc['F']['4_2'], percents_doc['F']['more'] ))

        #     # Quim: Writing an output file with the absolute number of differences dGexp - dGpred (Native)
        #     fan.write('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n' %( x,
        #                         absolute_nat['T']['1_4'], absolute_nat['T']['2_8'], absolute_nat['T']['4_2'], absolute_nat['T']['more'],
        #                         absolute_nat['R']['1_4'], absolute_nat['R']['2_8'], absolute_nat['R']['4_2'], absolute_nat['R']['more'],
        #                         absolute_nat['F']['1_4'], absolute_nat['F']['2_8'], absolute_nat['F']['4_2'], absolute_nat['F']['more'] ))

        #     # Quim: Writing an output file with the absolute number of differences dGexp - dGpred (All decoys)
        #     fad.write('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n' %( x,
        #                         absolute_doc['T']['1_4'], absolute_doc['T']['2_8'], absolute_doc['T']['4_2'], absolute_doc['T']['more'],
        #                         absolute_doc['R']['1_4'], absolute_doc['R']['2_8'], absolute_doc['R']['4_2'], absolute_doc['R']['more'],
        #                         absolute_doc['F']['1_4'], absolute_doc['F']['2_8'], absolute_doc['F']['4_2'], absolute_doc['F']['more'] ))


        #     distributions_table = 'distributions_{}.tsv'.format(x)
        #     fr = open(distributions_table, 'w')
        #     # <RIGID NATIVE>    <RIGID DOCK>  <FLEXIBLE NATIVE>   <FLEXIBLE DOCK>
        #     if len(R_predicted_dG) == len(F_predicted_dG) == len(R_predicted_dG_nat) == len(F_predicted_dG_nat):
        #         for num in xrange(len(R_predicted_dG)):
        #             fr.write('{}\t{}\t{}\t{}\n'.format( R_predicted_dG_nat[num], R_predicted_dG[num], F_predicted_dG_nat[num], F_predicted_dG[num] ))
        #     else:
        #         print(len(R_predicted_dG))
        #         print(len(F_predicted_dG))
        #         print(len(R_predicted_dG_nat))
        #         print(len(F_predicted_dG_nat))
        #     fr.close()


        # fpd.close()
        # fpn.close()
        # fad.close()
        # fan.close()

        R_docking_scores = [[i[enes_index['ES3DC']]] for i in R_docking]
        F_docking_scores = [[i[enes_index['ES3DC']]] for i in F_docking]
        print("Len R dG: {}  R scores: {}  F dG: {}  F scores: {}".format(len(R_dG), len(R_docking_scores), len(F_dG), len(F_docking_scores)))
        R_predicted_dG, R_experimental_dG, F_predicted_dG, F_experimental_dG, T_filt, R_filt, F_filt, _, _, corr = self.__traintest_linearmodel(R_docking_scores, F_docking_scores, R_dG, F_dG)


        # Print density plot (ES3DC) #
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
        plt.savefig(self.img_path+'/BA_ES3DC_Docking.eps', bbox_inches='tight', dpi=350)
        #plt.show()
        plt.clf()


        # Print scatter plot (ES3DC) #
        sns.set_context("paper", font_scale=1.5)
        sns.set_style("white")
        plt.ylim(-4,-19)
        plt.xlim(max(R_predicted_dG+F_predicted_dG), min(R_predicted_dG+F_predicted_dG))
        plt.xlabel(r'Predicted $\Delta$Gs (kcal mol$^{-1}$)')
        plt.ylabel(r'Experimental $\Delta$Gs (kcal mol$^{-1}$)')
        plt.scatter(np.array(R_filt['not_1_4']['pred']), np.array(R_filt['not_1_4']['exp']), c='#8e80fb', marker="o", edgecolor='none')
        plt.scatter(np.array(F_filt['not_1_4']['pred']), np.array(F_filt['not_1_4']['exp']), c='#ea968f', marker="v", edgecolor='none')
        plt.scatter(np.array(R_filt['1_4']['pred']), np.array(R_filt['1_4']['exp']), c='#1E02F7', marker="o", edgecolor='none')
        plt.scatter(np.array(F_filt['1_4']['pred']), np.array(F_filt['1_4']['exp']), c='#d62d20', marker="v", edgecolor='none')
        plt.savefig(self.img_path+'/BA_ES3DC_Docking_scatter_filt.eps', bbox_inches='tight', dpi=350)
        #plt.show()
        plt.clf()


        ####################################################################################

        # Print scatter plot of experimental dGs vs. docking scores (ES3DC) (with the regression line) #

        sns.set_context("paper", font_scale=1.5)
        sns.set_style("white")
        #plt.ylim(-4,-19)
        #plt.xlim(-8, -13.5)
        plt.xlabel(r'Docking scores')
        plt.ylabel(r'Experimental $\Delta$Gs (kcal mol$^{-1}$)')
        plt.scatter(np.array(R_docking_scores+F_docking_scores), np.array(R_dG+F_dG), c='blue', marker="o")
        
        X_plot = np.linspace(-8,18,100)
        plt.plot(X_plot, X_plot*corr[0] + corr[1], c='black')
        plt.plot(X_plot, X_plot*corr[0] + corr[1]+self.distance_reg_line, c='red')
        plt.plot(X_plot, X_plot*corr[0] + corr[1]-self.distance_reg_line, c='red')

        namefig = self.img_path+'/BA_ES3DC_Docking_scores_vs_dG_{}.eps'.format(self.distance_reg_line)
        plt.savefig(namefig, bbox_inches='tight', dpi=350)
        #plt.show()
        plt.clf()


        T_dock = R_docking_scores+F_docking_scores
        T_dG = R_dG+F_dG
        T_docking_scores = []
        for x in T_dock:
            T_docking_scores.append(x[0])
        zipped = zip(T_docking_scores, T_dG)
        num = 1
#        for score, dG in sorted(zipped, key=lambda tup: tup[0]):
#            print("{} --> SC: {}    DG: {}".format(num, score, dG))
#            num+=1

        # Take the docking values in R_dG and F_dG that are out of the area defined by the two lines of the best correlation

        R_docking_new = []
        F_docking_new = []
        R_dG_new = []
        F_dG_new = []
        zipped = zip(R_docking_scores, R_dG)
        for score, dG in zipped:
            score = score[0]
            upper_line = corr[0] * score + corr[1]+self.distance_reg_line
            lower_line = corr[0] * score + corr[1]-self.distance_reg_line
            if upper_line < dG or lower_line > dG:
                continue
            else: # Here we take only the points that are inside the area of the best correlation
                R_docking_new.append([score])
                R_dG_new.append(dG)
        zipped = zip(F_docking_scores, F_dG)
        for score, dG in zipped:
            score = score[0]
            upper_line = corr[0] * score + corr[1]+self.distance_reg_line
            lower_line = corr[0] * score + corr[1]-self.distance_reg_line
            if upper_line < dG or lower_line > dG:
                continue
            else: # Here we take only the points that are inside the area of the best correlation
                F_docking_new.append([score])
                F_dG_new.append(dG)


        # Print scatter plot of experimental dGs vs. docking scores removing the points (ES3DC) #
        sns.set_context("paper", font_scale=1.5)
        sns.set_style("white")
        #plt.ylim(-4,-19)
        #plt.xlim(-8, -13.5)
        plt.xlabel(r'Docking scores')
        plt.ylabel(r'Experimental $\Delta$Gs (kcal mol$^{-1}$)')
        plt.scatter(np.array(R_docking_new+F_docking_new), np.array(R_dG_new+F_dG_new), c='blue', marker="o")
        
        X_plot = np.linspace(-8,18,100)
        plt.plot(X_plot, X_plot*corr[0] + corr[1], c='black')
        plt.plot(X_plot, X_plot*corr[0] + corr[1]+self.distance_reg_line, c='red')
        plt.plot(X_plot, X_plot*corr[0] + corr[1]-self.distance_reg_line, c='red')

        namefig = self.img_path+'/BA_ES3DC_Docking_scores_vs_dG_points_removed_{}.eps'.format(self.distance_reg_line)
        plt.savefig(namefig, bbox_inches='tight', dpi=350)
        #plt.show()
        plt.clf()

        print("Len R new: {}  F new: {}".format(len(R_dG_new), len(F_dG_new)))
        # Calculating the predicted energies using the new pairs of values (pr = points removed)
        R_predicted_dG_pr, R_experimental_dG_pr, F_predicted_dG_pr, F_experimental_dG_pr, T_filt_pr, R_filt_pr, F_filt_pr, _, _, corr_pr = self.__traintest_linearmodel(R_docking_new, F_docking_new, R_dG_new, F_dG_new)

        # Print scatter plot with points removed (ES3DC) #
        sns.set_context("paper", font_scale=1.5)
        sns.set_style("white")
        plt.ylim(-4,-19)
        plt.xlim(max(R_predicted_dG+F_predicted_dG), min(R_predicted_dG+F_predicted_dG))
        plt.xlabel(r'Predicted $\Delta$Gs (kcal mol$^{-1}$)')
        plt.ylabel(r'Experimental $\Delta$Gs (kcal mol$^{-1}$)')
        plt.scatter(np.array(R_filt_pr['not_1_4']['pred']), np.array(R_filt_pr['not_1_4']['exp']), c='#8e80fb', marker="o", edgecolor='none')
        plt.scatter(np.array(F_filt_pr['not_1_4']['pred']), np.array(F_filt_pr['not_1_4']['exp']), c='#ea968f', marker="v", edgecolor='none')
        plt.scatter(np.array(R_filt_pr['1_4']['pred']), np.array(R_filt_pr['1_4']['exp']), c='#1E02F7', marker="o", edgecolor='none')
        plt.scatter(np.array(F_filt_pr['1_4']['pred']), np.array(F_filt_pr['1_4']['exp']), c='#d62d20', marker="v", edgecolor='none')
        namefig = self.img_path+'/BA_ES3DC_Docking_scatter_filt_points_removed_{}.eps'.format(self.distance_reg_line)
        plt.savefig(namefig, bbox_inches='tight', dpi=350)
        #plt.show()
        plt.clf()

        # Print density plot with points removed (ES3DC) #
        sns.set_context("paper", font_scale=1.5)
        sns.set_style("white")
        plt.ylim(-4,-19)
        plt.xlim(max(R_predicted_dG+F_predicted_dG), min(R_predicted_dG+F_predicted_dG))
        plt.xlabel(r'Predicted $\Delta$Gs (kcal mol$^{-1}$)')
        plt.ylabel(r'Experimental $\Delta$Gs (kcal mol$^{-1}$)')

        #sns.jointplot(np.array(R_predicted_dG+F_predicted_dG), np.array(R_experimental_dG+F_experimental_dG), kind="kde");
        sns.kdeplot(np.array(R_predicted_dG_pr), np.array(R_experimental_dG_pr), shade_lowest=False, cmap="Blues")
        sns.kdeplot(np.array(F_predicted_dG_pr), np.array(F_experimental_dG_pr), shade_lowest=False, cmap="Reds")
        #plt.scatter(R_predicted_dG, R_experimental_dG, color='b', marker='o', label='Rigid')
        #plt.scatter(F_predicted_dG, F_experimental_dG, color='r', marker='v', label='Flexible')
        #legend = plt.legend(loc='lower right', frameon=True)
        #legend.get_frame().set_facecolor('#E1E1E1')
        namefig = self.img_path+'/BA_ES3DC_Docking_density_points_removed_{}.eps'.format(self.distance_reg_line)
        plt.savefig(namefig, bbox_inches='tight', dpi=350)
        #plt.show()
        plt.clf()

        ####################################################################################



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
        
        # Quim: Defining some data structures
        Ty_diff_list = []
        Ry_diff_list = []
        Fy_diff_list = []
        T_filt = {'1_4': {'pred': [], 'exp': []},'2_8': {'pred': [], 'exp': []}, '4_2': {'pred': [], 'exp': []}, 'more': {'pred': [], 'exp': []}, 'not_1_4': {'pred': [], 'exp': []}}
        R_filt = {'1_4': {'pred': [], 'exp': []},'2_8': {'pred': [], 'exp': []}, '4_2': {'pred': [], 'exp': []}, 'more': {'pred': [], 'exp': []}, 'not_1_4': {'pred': [], 'exp': []}}
        F_filt = {'1_4': {'pred': [], 'exp': []},'2_8': {'pred': [], 'exp': []}, '4_2': {'pred': [], 'exp': []}, 'more': {'pred': [], 'exp': []}, 'not_1_4': {'pred': [], 'exp': []}}

        absolute = {'T': {'1_4': 0, '2_8': 0, '4_2': 0, 'more': 0}, 
            'R': {'1_4': 0, '2_8': 0, '4_2': 0, 'more': 0}, 
            'F': {'1_4': 0, '2_8': 0, '4_2': 0, 'more': 0}}


        for x in xrange(1000):
            # Quim notes: The scores would be the observations (X) / The dG would be the target values (y) http://scikit-learn.org/stable/tutorial/statistical_inference/supervised_learning.html
            # Quim notes: Now, we split the observations and known target values (experimental dG) in training and testing sets, and we build a model to predict dGs
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

            #print("Len R predict: {}  R test: {}  F predict: {}  F test: {}".format(len(Ry_predict), len(Ry_test), len(Fy_predict), len(Fy_test)))

            # Quim: Calculating the difference rigid predicted and rigid experimental values
            if len(Ry_predict) == len(Ry_test):
                for n in xrange(len(Ry_predict)):
                    Ry_diff = self.__calculate_diff_Gpred_Gexp(Ry_predict[n], Ry_test[n])
                    Ry_diff_list.append(Ry_diff)
                    if Ry_diff <= 1.4:
                        R_filt['1_4']['pred'].append(Ry_predict[n])
                        R_filt['1_4']['exp'].append(Ry_test[n])
                        absolute['R']['1_4'] += 1
                    else:
                        R_filt['not_1_4']['pred'].append(Ry_predict[n])
                        R_filt['not_1_4']['exp'].append(Ry_test[n])
                    if Ry_diff <= 2.8 and Ry_diff > 1.4:
                        R_filt['2_8']['pred'].append(Ry_predict[n])
                        R_filt['2_8']['exp'].append(Ry_test[n])
                        absolute['R']['2_8'] += 1
                    if Ry_diff <= 4.2  and Ry_diff > 2.8:
                        R_filt['4_2']['pred'].append(Ry_predict[n])
                        R_filt['4_2']['exp'].append(Ry_test[n])
                        absolute['R']['4_2'] += 1
                    if Ry_diff > 4.2:
                        R_filt['more']['pred'].append(Ry_predict[n])
                        R_filt['more']['exp'].append(Ry_test[n])
                        absolute['R']['more'] += 1
            else:
                print("ERROR: Ry predict, Rytest, Fy predict and Fy test not equal")
                print("Len R predict: {}  R test: {}".format(len(Ry_predict), len(Ry_test)))
                sys.exit(10)

            # Quim: Calculating the difference between the flexible predicted and the flexible experimental values
            if len(Fy_predict) == len(Fy_test):
                for n in xrange(len(Fy_predict)):
                    Fy_diff = self.__calculate_diff_Gpred_Gexp(Fy_predict[n], Fy_test[n])
                    Fy_diff_list.append(Fy_diff)
                    if Fy_diff <= 1.4:
                        F_filt['1_4']['pred'].append(Fy_predict[n])
                        F_filt['1_4']['exp'].append(Fy_test[n])
                        absolute['F']['1_4'] += 1
                    else:
                        F_filt['not_1_4']['pred'].append(Fy_predict[n])
                        F_filt['not_1_4']['exp'].append(Fy_test[n])
                    if Fy_diff <= 2.8 and Fy_diff > 1.4:
                        F_filt['2_8']['pred'].append(Fy_predict[n])
                        F_filt['2_8']['exp'].append(Fy_test[n])
                        absolute['F']['2_8'] += 1
                    if Fy_diff <= 4.2  and Fy_diff > 2.8:
                        F_filt['4_2']['pred'].append(Fy_predict[n])
                        F_filt['4_2']['exp'].append(Fy_test[n])
                        absolute['F']['4_2'] += 1
                    if Fy_diff > 4.2:
                        F_filt['more']['pred'].append(Fy_predict[n])
                        F_filt['more']['exp'].append(Fy_test[n])
                        absolute['F']['more'] += 1
            else:
                print("ERROR: Fy predict and Fy test are not equal")
                print("Len F predict: {}  F test: {}".format(len(Fy_predict), len(Fy_test)))
                sys.exit(10)

            # Quim: Calculating the difference between the total predicted and the total experimental values
            Ty_test = Ry_test+Fy_test
            if len(Ty_predict) == len(Ty_test):
                for n in xrange(len(Ty_predict)):
                    Ty_diff = self.__calculate_diff_Gpred_Gexp(Ty_predict[n], Ty_test[n])
                    Ty_diff_list.append(Ty_diff)
                    if Ty_diff <= 1.4:
                        T_filt['1_4']['pred'].append(Ty_predict[n])
                        T_filt['1_4']['exp'].append(Ty_test[n])
                        absolute['T']['1_4'] += 1
                    else:
                        T_filt['not_1_4']['pred'].append(Ty_predict[n])
                        T_filt['not_1_4']['exp'].append(Ty_test[n])
                    if Ty_diff <= 2.8 and Ty_diff > 1.4:
                        T_filt['2_8']['pred'].append(Ty_predict[n])
                        T_filt['2_8']['exp'].append(Ty_test[n])
                        absolute['T']['2_8'] += 1
                    if Ty_diff <= 4.2  and Ty_diff > 2.8:
                        T_filt['4_2']['pred'].append(Ty_predict[n])
                        T_filt['4_2']['exp'].append(Ty_test[n])
                        absolute['T']['4_2'] += 1
                    if Ty_diff > 4.2:
                        T_filt['more']['pred'].append(Ty_predict[n])
                        T_filt['more']['exp'].append(Ty_test[n])
                        absolute['T']['more'] += 1
            else:
                print("ERROR: Ty predict and Ty test are not equal")
                print("Len T predict: {}  T test: {}".format(len(Ty_predict), len(Ty_test)))
                sys.exit(10)


        # Quim: Calculating the percentages of differences of values below 1.4, 2.8 and 4.2
        percent_T_1_4 = float(len(T_filt['1_4']['pred']))/float(len(Ty_diff_list))*100
        percent_T_2_8 = float(len(T_filt['2_8']['pred']))/float(len(Ty_diff_list))*100
        percent_T_4_2 = float(len(T_filt['4_2']['pred']))/float(len(Ty_diff_list))*100
        percent_T_more = float(len(T_filt['more']['pred']))/float(len(Ty_diff_list))*100
        median_T = np.median(Ty_diff_list)
        print 'TOTAL -> Pearsonr: %.2f (std %.2f), RMSE: %.2f' % (np.nanmean(Tpearsons), np.nanstd(Tpearsons), np.nanmean(Tmses))
        print 'TOTAL -> Median absolute prediction error: %.2f' % (median_T)
        print 'TOTAL -> Prediction ratio 1.4: %.2f%%\t\tAbsolute: %d' % (percent_T_1_4, absolute['T']['1_4'])
        print 'TOTAL -> Prediction ratio 2.8: %.2f%%\t\tAbsolute: %d' % (percent_T_2_8, absolute['T']['2_8'])
        print 'TOTAL -> Prediction ratio 4.2: %.2f%%\t\tAbsolute: %d' % (percent_T_4_2, absolute['T']['4_2'])
        print 'TOTAL -> Prediction ratio > 4.2: %.2f%%\t\tAbsolute: %d' % (percent_T_more, absolute['T']['more'])
        print 'TOTAL -> Total predictions: %d' % (len(Ty_diff_list))

        percent_R_1_4 = float(len(R_filt['1_4']['pred']))/float(len(Ry_diff_list))*100
        percent_R_2_8 = float(len(R_filt['2_8']['pred']))/float(len(Ry_diff_list))*100
        percent_R_4_2 = float(len(R_filt['4_2']['pred']))/float(len(Ry_diff_list))*100
        percent_R_more = float(len(R_filt['more']['pred']))/float(len(Ry_diff_list))*100
        median_R = np.median(Ry_diff_list)
        print 'RIGID -> Pearsonr: %.2f (std %.2f), RMSE: %.2f' % (np.nanmean(Rpearsons), np.nanstd(Rpearsons), np.nanmean(Rmses))
        print 'RIGID -> Median absolute prediction error: %.2f' % (median_R)
        print 'RIGID -> Prediction ratio 1.4: %.2f%%' % (percent_R_1_4)
        print 'RIGID -> Prediction ratio 2.8: %.2f%%' % (percent_R_2_8)
        print 'RIGID -> Prediction ratio 4.2: %.2f%%' % (percent_R_4_2)
        print 'RIGID -> Prediction ratio > 4.2: %.2f%%' % (percent_R_more)
        print 'RIGID -> Total predictions: %d' % (len(Ry_diff_list))

        percent_F_1_4 = float(len(F_filt['1_4']['pred']))/float(len(Fy_diff_list))*100
        percent_F_2_8 = float(len(F_filt['2_8']['pred']))/float(len(Fy_diff_list))*100
        percent_F_4_2 = float(len(F_filt['4_2']['pred']))/float(len(Fy_diff_list))*100
        percent_F_more = float(len(F_filt['more']['pred']))/float(len(Fy_diff_list))*100
        median_F = np.median(Fy_diff_list)
        print 'FLEXI -> Pearsonr: %.2f (std %.2f), RMSE: %.2f' % (np.nanmean(Fpearsons), np.nanstd(Fpearsons), np.nanmean(Fmses))
        print 'FLEXI -> Median absolute prediction error: %.2f' % (median_F)
        print 'FLEXI -> Prediction ratio 1.4: %.2f%%' % (percent_F_1_4)
        print 'FLEXI -> Prediction ratio 2.8: %.2f%%' % (percent_F_2_8)
        print 'FLEXI -> Prediction ratio 4.2: %.2f%%' % (percent_F_4_2)
        print 'FLEXI -> Prediction ratio > 4.2: %.2f%%' % (percent_F_more)
        print 'FLEXI -> Total predictions: %d' % (len(Fy_diff_list))

        percents = {'T': {'1_4': percent_T_1_4, '2_8': percent_T_2_8, '4_2': percent_T_4_2, 'more': percent_T_more, 'median': median_T}, 
                    'R': {'1_4': percent_R_1_4, '2_8': percent_R_2_8, '4_2': percent_R_4_2, 'more': percent_R_more, 'median': median_R}, 
                    'F': {'1_4': percent_F_1_4, '2_8': percent_F_2_8, '4_2': percent_F_4_2, 'more': percent_F_more, 'median': median_F}}

        # Create model
        clf = linear_model.LinearRegression()
        clf.fit(R_scores+F_scores, R_dG+F_dG)
        print 'Coefs: %.10f*E, Intercept: %.3f' % (clf.coef_[0], clf.intercept_)
        corr_result = [clf.coef_[0], clf.intercept_]

        # print("NOW R SCORES")
        # print(R_scores)
        # print("NOW R_dG")
        # print(R_dG)

        return R_predicted_dG, R_experimental_dG, F_predicted_dG, F_experimental_dG, T_filt, R_filt, F_filt, percents, absolute, corr_result


    def __calculate_diff_Gpred_Gexp(self, G_pred, G_exp):
        '''
        Calculate the difference between dGpred and dGexp
        '''
        diff = abs(G_exp - G_pred)
        return diff


#--------------------------#
#          MAIN            #
#--------------------------#

if __name__ == "__main__":
   
    affinity_analysis = AffinityAnalysis() 
    affinity_analysis.run()
