import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
from scipy.stats import linregress
from scipy.stats.stats import pearsonr
import os, sys

def boxplot_ene(energies_file, energy_query):
    '''
    Plot the resultant energies in boxplots
    '''
    labels = ['FiberDock', 'aVdW', 'rVdW', 'ACE', 'inside', 'aElec', 'rElec', 'laElec', 'lrElec', 'hb',
              'EPair', 'ES3DC', 'ELOCAL', 'E3DC', 'E3D']
    energies = {}
    energies_fo = open(energies_file, 'r')
    for energy in energies_fo:
        energy = energy.strip().split('\t')[1:]
        for x in xrange(15):
            energies.setdefault(labels[x], {'N': [], 'I': [], 'P': [], 'E': []})
            energies[labels[x]]['N'].append(float(energy[x*4]))
            energies[labels[x]]['I'].append(float(energy[x*4+1]))
            energies[labels[x]]['P'].append(float(energy[x*4+2]))
            energies[labels[x]]['E'].append(float(energy[x*4+3]))

    N_mean = np.mean(energies[energy_query]['N'])
    I_mean = np.mean(energies[energy_query]['I'])
    P_mean = np.mean(energies[energy_query]['P'])
    E_mean = np.mean(energies[energy_query]['E'])
    plt.ylabel(energy_query)
    plt.xlabel('Decoy class')
    bp = plt.boxplot([energies[energy_query]['N'],
                 energies[energy_query]['I'],
                 energies[energy_query]['P'], 
                 energies[energy_query]['E']], 1, '')
    for box in bp['boxes']:
        box.set(color='k', linewidth=1)
    for box in bp['whiskers']:
        box.set(color='k', linewidth=1)
    for box in bp['caps']:
        box.set(color='k', linewidth=1)
    for box in bp['medians']:
        box.set(color='k', linewidth=1)
    plt.xticks([1, 2, 3, 4], ['Near-Native', 'Face-Face', 'Face-Back', 'Back-Back'])
    ax = plt.gca()
    if energy_query == 'E3D':
        text_mean = 'Average: NN %.0f, FF %.0f, FB %.0f, BB %.0f' % (N_mean, I_mean, P_mean, E_mean)
    else:
        text_mean = 'Average: NN %.2f, FF %.2f, FB %.2f, BB %.2f' % (N_mean, I_mean, P_mean, E_mean)
    ax.text(.012,.98, text_mean, ha='left', va='top', transform=ax.transAxes, size='smaller', bbox=dict(facecolor='#E1E1E1'))
    plt.savefig('figures/BP_%s.eps' % energy_query, bbox_inches='tight', dpi=350)
    plt.cla()

def scatter_ene(energies_file, energy_query):
    '''
    '''
    if energy_query == 'FiberDock':
        fmin = 0
    if energy_query == 'EPair':
        fmin = 40
    if energy_query == 'ES3DC':
        fmin = 44
    if energy_query == 'E3D':
        fmin = 56
    fmax = fmin+4

    ppi_energies = []
    ppi_errorene = []
    ppi_energies_fo = open(energies_file, 'r')
    for energy in ppi_energies_fo:
        energy = energy.strip().split('\t')
        mean_energy = [float(x) for x in energy[fmin+1:fmax+1]]
        if energy_query == 'FiberDock' and (mean_energy[0] > 500 or mean_energy[3] > 60):
            continue
        std_energy = [float(x) for x in energy[fmin+61:fmax+61]]
        ppi_energies.append(mean_energy)
        ppi_errorene.append(std_energy)
    ppi_energies_fo.close()
    ppi_energies = np.array(ppi_energies)
    ppi_errorene = np.array(ppi_errorene)

    # NN-FF
    pearson_r = pearsonr(ppi_energies[:,0], ppi_energies[:,1])
    slope, intercept, r_value, p_value, std_err = linregress(ppi_energies[:,0], ppi_energies[:,1])
    print 'Correlation NN-FF: slope %.2f, intercept %.2f' % (slope, intercept)
    plt.scatter(ppi_energies[:,0], ppi_energies[:,1], c='k', label='r=%.2f ; p-val=<0.01' % pearson_r[0])
    plt.errorbar(ppi_energies[:,0], ppi_energies[:,1], xerr=ppi_errorene[:,0], yerr=ppi_errorene[:,1], elinewidth=0.5, linestyle="None")
    plt.plot(ppi_energies[:,0], slope*ppi_energies[:,0]+intercept, c='k')
    plt.ylabel('Face-Face %s' % energy_query)
    plt.xlabel('Near-Native %s' % energy_query)
    legend = plt.legend(loc='upper left', frameon=True)
    legend.get_frame().set_facecolor('#E1E1E1')
    plt.savefig('figures/CorrNI_%s.eps' % energy_query, bbox_inches='tight', dpi=350)
    plt.cla()
    # NN-FB
    pearson_r = pearsonr(ppi_energies[:,0], ppi_energies[:,2])
    slope, intercept, r_value, p_value, std_err = linregress(ppi_energies[:,0], ppi_energies[:,2])
    print 'Correlation NN-FB: slope %.2f, intercept %.2f' % (slope, intercept)
    plt.scatter(ppi_energies[:,0], ppi_energies[:,2], c='k', label='r=%.2f ; p-val=<0.01' % pearson_r[0])
    plt.errorbar(ppi_energies[:,0], ppi_energies[:,2], xerr=ppi_errorene[:,0], yerr=ppi_errorene[:,2], elinewidth=0.5, linestyle="None")
    plt.plot(ppi_energies[:,0], slope*ppi_energies[:,0]+intercept, c='k')
    plt.ylabel('Face-Back %s' % energy_query)
    plt.xlabel('Near-Native %s' % energy_query)
    legend = plt.legend(loc='upper left', frameon=True)
    legend.get_frame().set_facecolor('#E1E1E1')
    plt.savefig('figures/CorrNP_%s.eps' % energy_query, bbox_inches='tight', dpi=350)
    plt.cla()
    # NN-BB
    pearson_r = pearsonr(ppi_energies[:,0], ppi_energies[:,3])
    slope, intercept, r_value, p_value, std_err = linregress(ppi_energies[:,0], ppi_energies[:,3])
    print 'Correlation NN-BB: slope %.2f, intercept %.2f' % (slope, intercept)
    plt.scatter(ppi_energies[:,0], ppi_energies[:,3], c='k', label='r=%.2f ; p-val=<0.01' % pearson_r[0])
    plt.errorbar(ppi_energies[:,0], ppi_energies[:,3], xerr=ppi_errorene[:,0], yerr=ppi_errorene[:,3], elinewidth=0.5, linestyle="None")
    plt.plot(ppi_energies[:,0], slope*ppi_energies[:,0]+intercept, c='k')
    plt.ylabel('Back-Back %s' % energy_query)
    plt.xlabel('Near-Native %s' % energy_query)
    legend = plt.legend(loc='upper left', frameon=True)
    legend.get_frame().set_facecolor('#E1E1E1')
    plt.savefig('figures/CorrNE_%s.eps' % energy_query, bbox_inches='tight', dpi=350)
    plt.cla()
    # FF-FB
    slope, intercept, r_value, p_value, std_err = linregress(ppi_energies[:,1], ppi_energies[:,2])
    print 'Correlation FF-FB: slope %.2f, intercept %.2f' % (slope, intercept)
    pearson_r = pearsonr(ppi_energies[:,1], ppi_energies[:,2])
    plt.scatter(ppi_energies[:,1], ppi_energies[:,2], c='k', label='r=%.2f ; p-val=<0.01' % pearson_r[0])
    plt.errorbar(ppi_energies[:,1], ppi_energies[:,2], xerr=ppi_errorene[:,1], yerr=ppi_errorene[:,2], elinewidth=0.5, linestyle="None")
    plt.plot(ppi_energies[:,1], slope*ppi_energies[:,1]+intercept, c='k')
    plt.ylabel('Face-Back %s' % energy_query)
    plt.xlabel('Face-Face %s' % energy_query)
    legend = plt.legend(loc='upper left', frameon=True)
    legend.get_frame().set_facecolor('#E1E1E1')
    plt.savefig('figures/CorrIP_%s.eps' % energy_query, bbox_inches='tight', dpi=350)
    plt.cla()
    # FF-BB
    slope, intercept, r_value, p_value, std_err = linregress(ppi_energies[:,1], ppi_energies[:,3])
    print 'Correlation FF-BB: slope %.2f, intercept %.2f' % (slope, intercept)
    pearson_r = pearsonr(ppi_energies[:,1], ppi_energies[:,3])
    plt.scatter(ppi_energies[:,1], ppi_energies[:,3], c='k', label='r=%.2f ; p-val=<0.01' % pearson_r[0])
    plt.errorbar(ppi_energies[:,1], ppi_energies[:,3], xerr=ppi_errorene[:,1], yerr=ppi_errorene[:,3], elinewidth=0.5, linestyle="None")
    plt.plot(ppi_energies[:,1], slope*ppi_energies[:,1]+intercept, c='k')
    plt.ylabel('Back-Back %s' % energy_query)
    plt.xlabel('Face-Face %s' % energy_query)
    legend = plt.legend(loc='upper left', frameon=True)
    legend.get_frame().set_facecolor('#E1E1E1')
    plt.savefig('figures/CorrIE_%s.eps' % energy_query, bbox_inches='tight', dpi=350)
    plt.cla()
    # FB-BB
    slope, intercept, r_value, p_value, std_err = linregress(ppi_energies[:,2], ppi_energies[:,3])
    print 'Correlation FB-BB: slope %.2f, intercept %.2f' % (slope, intercept)
    pearson_r = pearsonr(ppi_energies[:,2], ppi_energies[:,3])
    plt.scatter(ppi_energies[:,2], ppi_energies[:,3], c='k', label='r=%.2f ; p-val=<0.01' % pearson_r[0])
    plt.errorbar(ppi_energies[:,2], ppi_energies[:,3], xerr=ppi_errorene[:,2], yerr=ppi_errorene[:,3], elinewidth=0.5, linestyle="None")
    plt.plot(ppi_energies[:,2], slope*ppi_energies[:,2]+intercept, c='k')
    plt.ylabel('Back-Back %s' % energy_query)
    plt.xlabel('Face-Back %s' % energy_query)
    legend = plt.legend(loc='upper left', frameon=True)
    legend.get_frame().set_facecolor('#E1E1E1')
    plt.savefig('figures/CorrPE_%s.eps' % energy_query, bbox_inches='tight', dpi=350)
    plt.cla()

def get_affinities():
    '''
    Get the complex affinities of Affinity Benchmark 2
    '''
    affinities = {}
    affbench_path = os.path.expanduser('~/Databases/AffinityBenchmark')
    ba_file = os.path.join(affbench_path, 'affinity_benchmark_all_dG.tsv')
    ba_fo = open(ba_file, 'r')
    for ppi in ba_fo:
        ppi = ppi.strip('\n').split('\t')
        ppi_name = ppi[0]
        dG = float(ppi[3])
        affinities[ppi_name] = dG
    ba_fo.close()
    return affinities

def affinities_corr(energies_file, energy_query):
    '''
    '''
    affinities = get_affinities()
    if energy_query == 'FiberDock':
        fmin = 0
    if energy_query == 'EPair':
        fmin = 40
    if energy_query == 'ES3DC':
        fmin = 44
    if energy_query == 'E3D':
        fmin = 56
    fmax = fmin+4

    ppi_affinities = []
    ppi_energies = []
    ppi_energies_fo = open(energies_file, 'r')
    for energy in ppi_energies_fo:
        ppi_name = energy[0:4]
        energy = [float(x) for x in energy.strip().split('\t')[fmin+1:fmax+1]]
        if energy_query == 'FiberDock' and (energy[0] > 500 or energy[3] > 60):
            continue
        if ppi_name in affinities:
            ppi_energies.append(energy)
            ppi_affinities.append(affinities[ppi_name])
    ppi_energies_fo.close()
    ppi_energies = np.array(ppi_energies)
    ppi_affinities = np.array(ppi_affinities)

    print 'Affinity NN: %.3f %.3e' % pearsonr(ppi_energies[:,0], ppi_affinities)
    print 'Affinity FF: %.3f %.3e' % pearsonr(ppi_energies[:,1], ppi_affinities)
    print 'Affinity FB: %.3f %.3e' % pearsonr(ppi_energies[:,2], ppi_affinities)
    print 'Affinity BB: %.3f %.3e' % pearsonr(ppi_energies[:,3], ppi_affinities)
    
#--------------------------#
#          MAIN            #
#--------------------------#

if __name__ == "__main__":

    enes_locations = os.path.abspath('scores_locations_E+OX_std.txt')
    sns.set_context("paper", font_scale=1.5)
    sns.set_style("white")
    for ene_type in ['FiberDock', 'E3D', 'EPair', 'ES3DC']:
        print '---------------- '+ene_type+' ----------------'
        boxplot_ene(enes_locations, ene_type)
        scatter_ene(enes_locations, ene_type)
        montage_crr = 'montage -tile 2x3 -geometry +0+0 -border 2 '
        montage_crr += 'figures/CorrNI_%s.eps figures/CorrIP_%s.eps ' % (ene_type, ene_type)
        montage_crr += 'figures/CorrNP_%s.eps figures/CorrPE_%s.eps ' % (ene_type, ene_type)
        montage_crr += 'figures/CorrNE_%s.eps figures/CorrIE_%s.eps ' % (ene_type, ene_type)
        montage_crr += 'figures/Scatters_%s.png' % ene_type
        os.system(montage_crr)
        affinities_corr(enes_locations, ene_type)
