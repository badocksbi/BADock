import cPickle
import os
from os.path import join
import re


def main():
    """
    Performs several analyses related with the calculation of the percentages for the affinity analysis
    """
    results_dir = "results_quim" # Introduce your results folder path
    try:
        os.stat("temp")
    except:
        os.mkdir("temp")

    # For TOTAL
    dump_pred = join(results_dir,"predictions.pcl")
    dump_dGexp = join(results_dir,"dGexp.pcl")
    dump_pred_CCharPPI = join(results_dir,"predictions_CCharPPI.pcl")
    dump_dGexp_CCharPPI = join(results_dir,"dGexp_CCharPPI.pcl")

    total_predictions = cPickle.load(open(dump_pred))
    total_dGs = cPickle.load(open(dump_dGexp))
    total_predictions_CCharPPI = cPickle.load(open(dump_pred_CCharPPI))
    total_dGs_CCharPPI = cPickle.load(open(dump_dGexp_CCharPPI))

    # For RIGID
    dump_pred_R = join(results_dir,"R_predictions.pcl")
    dump_dGexp_R = join(results_dir,"R_dGexp.pcl")
    dump_pred_CCharPPI_R = join(results_dir,"R_predictions_CCharPPI.pcl")
    dump_dGexp_CCharPPI_R = join(results_dir,"R_dGexp_CCharPPI.pcl")

    R_predictions = cPickle.load(open(dump_pred_R))
    R_dG = cPickle.load(open(dump_dGexp_R))
    R_predictions_CCharPPI = cPickle.load(open(dump_pred_CCharPPI_R))
    R_dG_CCharPPI = cPickle.load(open(dump_dGexp_CCharPPI_R))

    # For FLEXIBLE
    dump_pred_F = join(results_dir,"F_predictions.pcl")
    dump_dGexp_F = join(results_dir,"F_dGexp.pcl")
    dump_pred_CCharPPI_F = join(results_dir,"F_predictions_CCharPPI.pcl")
    dump_dGexp_CCharPPI_F = join(results_dir,"F_dGexp_CCharPPI.pcl")

    F_predictions = cPickle.load(open(dump_pred_F))
    F_dG = cPickle.load(open(dump_dGexp_F))
    F_predictions_CCharPPI = cPickle.load(open(dump_pred_CCharPPI_F))
    F_dG_CCharPPI = cPickle.load(open(dump_dGexp_CCharPPI_F))


    print('\n############# TOTAL ANALYSIS #############')
    T_filtered, T_percents = calculate_percents_2_rows(total_dGs, total_predictions)
    T_filtered_CCharPPI, T_percents_CCharPPI = calculate_percents_2_rows(total_dGs_CCharPPI, total_predictions_CCharPPI, CCharPPI=True)

    print('\n############# RIGID ANALYSIS #############')
    R_filtered, R_percents = calculate_percents_2_rows(R_dG, R_predictions)
    R_filtered_CCharPPI, R_percents_CCharPPI = calculate_percents_2_rows(R_dG_CCharPPI, R_predictions_CCharPPI, CCharPPI=True)

    print('\n############# FLEXIBLE ANALYSIS #############')
    F_filtered, F_percents = calculate_percents_2_rows(F_dG, F_predictions)
    F_filtered_CCharPPI, F_percents_CCharPPI = calculate_percents_2_rows(F_dG_CCharPPI, F_predictions_CCharPPI, CCharPPI=True)

    # Writing table of results
    output_file = join(results_dir,"percents_individual_predictions_complete_2_rows.tsv")
    outfile = open(output_file, 'w')
    print_table_2_rows(T_percents, T_filtered, T_percents_CCharPPI, T_filtered_CCharPPI, R_percents, R_filtered, R_percents_CCharPPI, R_filtered_CCharPPI, F_percents, F_filtered, F_percents_CCharPPI, F_filtered_CCharPPI, outfile)
    outfile.close()

    # Perform Mann Whitney
    mw_results_T = calculate_mann_whitney_2_rows(T_filtered)
    mw_results_R = calculate_mann_whitney_2_rows(R_filtered)
    mw_results_F = calculate_mann_whitney_2_rows(F_filtered)
    print mw_results_T

    # Writing table of Mann Whitney results
    mann_whitney_file = join(results_dir,"mann_whitney_table_2_rows.tsv")
    mwfile = open(mann_whitney_file, 'w')
    print_mann_whitney_2_rows(mw_results_T, mw_results_R, mw_results_F, mwfile)
    mwfile.close()

    return


def calculate_percents_cumulative(dGs_experimental, predictions, CCharPPI=False):
    '''
    Calculates the percentages
    '''

    filtered = {}
    percents = {}

    for potential in predictions:

        print '\n==== %s ====' %(potential.upper())
        filtered.setdefault(potential, {})
        percents.setdefault(potential, {})

        for type_struc in predictions[potential]:

            if CCharPPI == True:
                if type_struc == 'docking':
                    continue

            filtered[potential].setdefault(type_struc, { '1_4': {'pred': [], 'exp': []},'2_8': {'pred': [], 'exp': []}, '4_2': {'pred': [], 'exp': []}, 'more': {'pred': [], 'exp': []} })
            percents[potential].setdefault(type_struc, { '1_4': 0, '2_8': 0, '4_2': 0, 'more': 0 })
            zipped = zip(dGs_experimental, predictions[potential][type_struc])

            for dGexp, pred in zipped:

                diff = calculate_diff_Gexp_Gpred(dGexp, pred)

                if diff <= 1.4:
                    filtered[potential][type_struc]['1_4']['pred'].append(pred)
                    filtered[potential][type_struc]['1_4']['exp'].append(dGexp)
                if diff <= 2.8:
                    filtered[potential][type_struc]['2_8']['pred'].append(pred)
                    filtered[potential][type_struc]['2_8']['exp'].append(dGexp)
                if diff <= 4.2:
                    filtered[potential][type_struc]['4_2']['pred'].append(pred)
                    filtered[potential][type_struc]['4_2']['exp'].append(dGexp)
                if diff > 4.2:
                    filtered[potential][type_struc]['more']['pred'].append(pred)
                    filtered[potential][type_struc]['more']['exp'].append(dGexp)

            percent_1_4 = float(len(filtered[potential][type_struc]['1_4']['pred'])) / float(len(predictions[potential][type_struc])) * 100
            percent_2_8 = float(len(filtered[potential][type_struc]['2_8']['pred'])) / float(len(predictions[potential][type_struc])) * 100
            percent_4_2 = float(len(filtered[potential][type_struc]['4_2']['pred'])) / float(len(predictions[potential][type_struc])) * 100
            percent_more = float(len(filtered[potential][type_struc]['more']['pred'])) / float(len(predictions[potential][type_struc])) * 100
            percents[potential][type_struc]['1_4'] = percent_1_4
            percents[potential][type_struc]['2_8'] = percent_2_8
            percents[potential][type_struc]['4_2'] = percent_4_2
            percents[potential][type_struc]['more'] = percent_more
            print '---> %s <---' %(type_struc)
            print 'Prediction ratio 1.4:\t%.2f%%\tAbsolute: %d' % (percent_1_4, len(filtered[potential][type_struc]['1_4']['pred']))
            print 'Prediction ratio 2.8:\t%.2f%%\tAbsolute: %d' % (percent_2_8, len(filtered[potential][type_struc]['2_8']['pred']))
            print 'Prediction ratio 4.2:\t%.2f%%\tAbsolute: %d' % (percent_4_2, len(filtered[potential][type_struc]['4_2']['pred']))
            print 'Prediction ratio > 4.2:\t%.2f%%\tAbsolute: %d' % (percent_more, len(filtered[potential][type_struc]['more']['pred']))

    return filtered, percents


def calculate_percents_non_cumulative(dGs_experimental, predictions, CCharPPI=False):
    '''
    Calculates the percentages
    '''

    filtered = {}
    percents = {}

    for potential in predictions:

        print '\n==== %s ====' %(potential.upper())
        filtered.setdefault(potential, {})
        percents.setdefault(potential, {})

        for type_struc in predictions[potential]:

            if CCharPPI == True:
                if type_struc == 'docking':
                    continue

            filtered[potential].setdefault(type_struc, { '1_4': {'pred': [], 'exp': []},'2_8': {'pred': [], 'exp': []}, '4_2': {'pred': [], 'exp': []}, 'more': {'pred': [], 'exp': []} })
            percents[potential].setdefault(type_struc, { '1_4': 0, '2_8': 0, '4_2': 0, 'more': 0 })
            zipped = zip(dGs_experimental, predictions[potential][type_struc])

            for dGexp, pred in zipped:

                diff = calculate_diff_Gexp_Gpred(dGexp, pred)

                if diff <= 1.4:
                    filtered[potential][type_struc]['1_4']['pred'].append(pred)
                    filtered[potential][type_struc]['1_4']['exp'].append(dGexp)
                if diff <= 2.8 and diff > 1.4:
                    filtered[potential][type_struc]['2_8']['pred'].append(pred)
                    filtered[potential][type_struc]['2_8']['exp'].append(dGexp)
                if diff <= 4.2  and diff > 2.8:
                    filtered[potential][type_struc]['4_2']['pred'].append(pred)
                    filtered[potential][type_struc]['4_2']['exp'].append(dGexp)
                if diff > 4.2:
                    filtered[potential][type_struc]['more']['pred'].append(pred)
                    filtered[potential][type_struc]['more']['exp'].append(dGexp)

            percent_1_4 = float(len(filtered[potential][type_struc]['1_4']['pred'])) / float(len(predictions[potential][type_struc])) * 100
            percent_2_8 = float(len(filtered[potential][type_struc]['2_8']['pred'])) / float(len(predictions[potential][type_struc])) * 100
            percent_4_2 = float(len(filtered[potential][type_struc]['4_2']['pred'])) / float(len(predictions[potential][type_struc])) * 100
            percent_more = float(len(filtered[potential][type_struc]['more']['pred'])) / float(len(predictions[potential][type_struc])) * 100
            percents[potential][type_struc]['1_4'] = percent_1_4
            percents[potential][type_struc]['2_8'] = percent_2_8
            percents[potential][type_struc]['4_2'] = percent_4_2
            percents[potential][type_struc]['more'] = percent_more
            print '---> %s <---' %(type_struc)
            print 'Prediction ratio 1.4:\t%.2f%%\tAbsolute: %d' % (percent_1_4, len(filtered[potential][type_struc]['1_4']['pred']))
            print 'Prediction ratio 2.8:\t%.2f%%\tAbsolute: %d' % (percent_2_8, len(filtered[potential][type_struc]['2_8']['pred']))
            print 'Prediction ratio 4.2:\t%.2f%%\tAbsolute: %d' % (percent_4_2, len(filtered[potential][type_struc]['4_2']['pred']))
            print 'Prediction ratio > 4.2:\t%.2f%%\tAbsolute: %d' % (percent_more, len(filtered[potential][type_struc]['more']['pred']))

    return filtered, percents


def calculate_percents_2_rows(dGs_experimental, predictions, CCharPPI=False):
    '''
    Calculates the percentages
    '''

    filtered = {}
    percents = {}

    for potential in predictions:

        print '\n==== %s ====' %(potential.upper())
        filtered.setdefault(potential, {})
        percents.setdefault(potential, {})

        for type_struc in predictions[potential]:

            if CCharPPI == True:
                if type_struc == 'docking':
                    continue

            filtered[potential].setdefault(type_struc, { '2_8': {'pred': [], 'exp': []}, 'more': {'pred': [], 'exp': []} })
            percents[potential].setdefault(type_struc, { '2_8': 0, 'more': 0 })
            zipped = zip(dGs_experimental, predictions[potential][type_struc])

            for dGexp, pred in zipped:

                diff = calculate_diff_Gexp_Gpred(dGexp, pred)

                if diff <= 2.8:
                    filtered[potential][type_struc]['2_8']['pred'].append(pred)
                    filtered[potential][type_struc]['2_8']['exp'].append(dGexp)
                if diff > 2.8:
                    filtered[potential][type_struc]['more']['pred'].append(pred)
                    filtered[potential][type_struc]['more']['exp'].append(dGexp)

            percent_2_8 = float(len(filtered[potential][type_struc]['2_8']['pred'])) / float(len(predictions[potential][type_struc])) * 100
            percent_more = float(len(filtered[potential][type_struc]['more']['pred'])) / float(len(predictions[potential][type_struc])) * 100
            percents[potential][type_struc]['2_8'] = percent_2_8
            percents[potential][type_struc]['more'] = percent_more
            print '---> %s <---' %(type_struc)
            print 'Prediction ratio 2.8:\t%.2f%%\tAbsolute: %d' % (percent_2_8, len(filtered[potential][type_struc]['2_8']['pred']))
            print 'Prediction ratio > 2.8:\t%.2f%%\tAbsolute: %d' % (percent_more, len(filtered[potential][type_struc]['more']['pred']))

    return filtered, percents


def calculate_diff_Gexp_Gpred(G_exp, G_pred):
    '''
    Calculate the difference between dGexp and dGpred
    '''
    diff = abs(G_exp - G_pred)
    return diff


def print_table(percents, filtered, percents_CCharPPI, filtered_CCharPPI, percents_R, filtered_R, percents_R_CCharPPI, filtered_R_CCharPPI, percents_F, filtered_F, percents_F_CCharPPI, filtered_F_CCharPPI, outfile):
    '''
    Prints the table with the results
    '''

    for potential in ['FiberDock', 'aVdW', 'rVdW', 'aElec', 'rElec', 'laElec', 'lrElec', 'HB', 'EPAIR', 'ES3DC', 'E3D']:

        potential_name = potential
        if potential == 'ROSSETADOCK':
            potential_name = 'RosettaDock'

        for percent in ['1.4', '2.8', '4.2', '>4.2']:
            if percent == '1.4':
                outfile.write("{}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\n".format( potential_name, '<=1.4', percents[potential]['native']['1_4'], len(filtered[potential]['native']['1_4']['pred']), percents[potential]['docking']['1_4'], len(filtered[potential]['docking']['1_4']['pred']),
                                                                                                                        percents_R[potential]['native']['1_4'], len(filtered_R[potential]['native']['1_4']['pred']), percents_R[potential]['docking']['1_4'], len(filtered_R[potential]['docking']['1_4']['pred']),
                                                                                                                        percents_F[potential]['native']['1_4'], len(filtered_F[potential]['native']['1_4']['pred']), percents_F[potential]['docking']['1_4'], len(filtered_F[potential]['docking']['1_4']['pred']) ))
            elif percent == '2.8':
                outfile.write("{}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\n".format( '', '<=2.8', percents[potential]['native']['2_8'], len(filtered[potential]['native']['2_8']['pred']), percents[potential]['docking']['2_8'], len(filtered[potential]['docking']['2_8']['pred']),
                                                                                                                        percents_R[potential]['native']['2_8'], len(filtered_R[potential]['native']['2_8']['pred']), percents_R[potential]['docking']['2_8'], len(filtered_R[potential]['docking']['2_8']['pred']),
                                                                                                                        percents_F[potential]['native']['2_8'], len(filtered_F[potential]['native']['2_8']['pred']), percents_F[potential]['docking']['2_8'], len(filtered_F[potential]['docking']['2_8']['pred']) ))
            elif percent == '4.2':
                outfile.write("{}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\n".format( '', '<=4.2', percents[potential]['native']['4_2'], len(filtered[potential]['native']['4_2']['pred']), percents[potential]['docking']['4_2'], len(filtered[potential]['docking']['4_2']['pred']),
                                                                                                                        percents_R[potential]['native']['4_2'], len(filtered_R[potential]['native']['4_2']['pred']), percents_R[potential]['docking']['4_2'], len(filtered_R[potential]['docking']['4_2']['pred']),
                                                                                                                        percents_F[potential]['native']['4_2'], len(filtered_F[potential]['native']['4_2']['pred']), percents_F[potential]['docking']['4_2'], len(filtered_F[potential]['docking']['4_2']['pred']) ))
            elif percent == '>4.2':
                outfile.write("{}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\n".format( '', '>4.2', percents[potential]['native']['more'], len(filtered[potential]['native']['more']['pred']), percents[potential]['docking']['more'], len(filtered[potential]['docking']['more']['pred']),
                                                                                                                        percents_R[potential]['native']['more'], len(filtered_R[potential]['native']['more']['pred']), percents_R[potential]['docking']['more'], len(filtered_R[potential]['docking']['more']['pred']),
                                                                                                                        percents_F[potential]['native']['more'], len(filtered_F[potential]['native']['more']['pred']), percents_F[potential]['docking']['more'], len(filtered_F[potential]['docking']['more']['pred']) ))

    for potential in ['ZRANK', 'ZRANK2', 'ROSSETADOCK', 'PYDOCK', 'PISA', 'PIE', 'SIPPER']:

        potential_name = potential
        if potential == 'PYDOCK':
            potential_name = 'PyDock'

        for percent in ['1.4', '2.8', '4.2', '>4.2']:
            if percent == '1.4':
                outfile.write("{}\t{}\t{:.2f}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\n".format( potential_name, '<=1.4', percents_CCharPPI[potential]['native']['1_4'], len(filtered_CCharPPI[potential]['native']['1_4']['pred']), '-', '-',
                                                                                                                        percents_R_CCharPPI[potential]['native']['1_4'], len(filtered_R_CCharPPI[potential]['native']['1_4']['pred']), '-', '-',
                                                                                                                        percents_F_CCharPPI[potential]['native']['1_4'], len(filtered_F_CCharPPI[potential]['native']['1_4']['pred']), '-', '-' ))
            elif percent == '2.8':
                outfile.write("{}\t{}\t{:.2f}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\n".format( '', '<=2.8', percents_CCharPPI[potential]['native']['2_8'], len(filtered_CCharPPI[potential]['native']['2_8']['pred']), '-', '-',
                                                                                                                        percents_R_CCharPPI[potential]['native']['2_8'], len(filtered_R_CCharPPI[potential]['native']['2_8']['pred']), '-', '-',
                                                                                                                        percents_F_CCharPPI[potential]['native']['2_8'], len(filtered_F_CCharPPI[potential]['native']['2_8']['pred']), '-', '-' ))
            elif percent == '4.2':
                outfile.write("{}\t{}\t{:.2f}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\n".format( '', '<=4.2', percents_CCharPPI[potential]['native']['4_2'], len(filtered_CCharPPI[potential]['native']['4_2']['pred']), '-', '-',
                                                                                                                        percents_R_CCharPPI[potential]['native']['4_2'], len(filtered_R_CCharPPI[potential]['native']['4_2']['pred']), '-', '-',
                                                                                                                        percents_F_CCharPPI[potential]['native']['4_2'], len(filtered_F_CCharPPI[potential]['native']['4_2']['pred']), '-', '-' ))
            elif percent == '>4.2':
                outfile.write("{}\t{}\t{:.2f}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\n".format( '', '>4.2', percents_CCharPPI[potential]['native']['more'], len(filtered_CCharPPI[potential]['native']['more']['pred']), '-', '-',
                                                                                                                        percents_R_CCharPPI[potential]['native']['more'], len(filtered_R_CCharPPI[potential]['native']['more']['pred']), '-', '-',
                                                                                                                        percents_F_CCharPPI[potential]['native']['more'], len(filtered_F_CCharPPI[potential]['native']['more']['pred']), '-', '-' ))

    return


def print_table_2_rows(percents, filtered, percents_CCharPPI, filtered_CCharPPI, percents_R, filtered_R, percents_R_CCharPPI, filtered_R_CCharPPI, percents_F, filtered_F, percents_F_CCharPPI, filtered_F_CCharPPI, outfile):
    '''
    Prints the table with the results
    '''

    for potential in ['FiberDock', 'aVdW', 'rVdW', 'aElec', 'rElec', 'laElec', 'lrElec', 'HB', 'EPAIR', 'ES3DC', 'E3D']:

        potential_name = potential
        if potential == 'ROSSETADOCK':
            potential_name = 'RosettaDock'

        for percent in ['2.8', '>2.8']:
            if percent == '2.8':
                outfile.write("{}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\n".format( potential_name, '<=2.8', percents[potential]['native']['2_8'], len(filtered[potential]['native']['2_8']['pred']), percents[potential]['docking']['2_8'], len(filtered[potential]['docking']['2_8']['pred']),
                                                                                                                        percents_R[potential]['native']['2_8'], len(filtered_R[potential]['native']['2_8']['pred']), percents_R[potential]['docking']['2_8'], len(filtered_R[potential]['docking']['2_8']['pred']),
                                                                                                                        percents_F[potential]['native']['2_8'], len(filtered_F[potential]['native']['2_8']['pred']), percents_F[potential]['docking']['2_8'], len(filtered_F[potential]['docking']['2_8']['pred']) ))
            elif percent == '>2.8':
                outfile.write("{}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\n".format( '', '>2.8', percents[potential]['native']['more'], len(filtered[potential]['native']['more']['pred']), percents[potential]['docking']['more'], len(filtered[potential]['docking']['more']['pred']),
                                                                                                                        percents_R[potential]['native']['more'], len(filtered_R[potential]['native']['more']['pred']), percents_R[potential]['docking']['more'], len(filtered_R[potential]['docking']['more']['pred']),
                                                                                                                        percents_F[potential]['native']['more'], len(filtered_F[potential]['native']['more']['pred']), percents_F[potential]['docking']['more'], len(filtered_F[potential]['docking']['more']['pred']) ))

    for potential in ['ZRANK', 'ZRANK2', 'ROSSETADOCK', 'PYDOCK', 'PISA', 'PIE', 'SIPPER']:

        potential_name = potential
        if potential == 'PYDOCK':
            potential_name = 'PyDock'

        for percent in ['2.8', '>2.8']:
            if percent == '2.8':
                outfile.write("{}\t{}\t{:.2f}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\n".format( potential_name, '<=2.8', percents_CCharPPI[potential]['native']['2_8'], len(filtered_CCharPPI[potential]['native']['2_8']['pred']), '-', '-',
                                                                                                                        percents_R_CCharPPI[potential]['native']['2_8'], len(filtered_R_CCharPPI[potential]['native']['2_8']['pred']), '-', '-',
                                                                                                                        percents_F_CCharPPI[potential]['native']['2_8'], len(filtered_F_CCharPPI[potential]['native']['2_8']['pred']), '-', '-' ))
            elif percent == '>2.8':
                outfile.write("{}\t{}\t{:.2f}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\n".format( '', '>2.8', percents_CCharPPI[potential]['native']['more'], len(filtered_CCharPPI[potential]['native']['more']['pred']), '-', '-',
                                                                                                                        percents_R_CCharPPI[potential]['native']['more'], len(filtered_R_CCharPPI[potential]['native']['more']['pred']), '-', '-',
                                                                                                                        percents_F_CCharPPI[potential]['native']['more'], len(filtered_F_CCharPPI[potential]['native']['more']['pred']), '-', '-' ))

    return


def calculate_mann_whitney(filtered):
    '''
    Calculates the Mann Whitney test for the pairs of Native - All decoys values
    '''

    mw_results = {}

    for potential in filtered:

        mw_results.setdefault(potential, {})

        for limit in ('1_4', '2_8', '4_2', 'more'):

            mw_results[potential].setdefault(limit, {})

            native_list = filtered[potential]['native'][limit]['pred']
            docking_list = filtered[potential]['docking'][limit]['pred']
            temp_file = 'temp/temp.R'
            temp = open(temp_file, 'w')

            for values_list, type_list in ((native_list, 'native'), (docking_list, 'docking')):
                temp.write('{} <- c('.format(type_list))
                for x in xrange(len(values_list)):
                    if x == len(values_list) - 1:
                        temp.write('{})\n'.format(values_list[x]))
                    else:
                        temp.write('{}, '.format(values_list[x]))
            temp.write('mwt <- wilcox.test(native, docking, conf.int = TRUE)\n')
            temp.write('mwt$p.value\n')
            temp.write('mwt$estimate\n')
            temp.write('mwt$conf.int[1:2]\n')
            temp.close()

            command = '/soft/R/R-3.3.0/bin/Rscript temp/temp.R > temp/temp.out'
            os.system(command)

            out = open('temp/temp.out', 'r')
            num_line = 1
            for line in out:
                if num_line == 1:
                    line1 = line.strip()
                if num_line == 3:
                    line3 = line.strip()
                if num_line == 4:
                    line4 = line.strip()
                num_line+=1
            pval = line1.split('[1] ')[1]
            mw_results[potential][limit]["pvalue"] = pval
            print("PVAL: {}".format(pval))
            estimate = ""
            estimate_pattern = re.compile('\s*([-+]?[0-9]*\.?[0-9]*)')
            m = estimate_pattern.search(line3)
            if m:
                estimate = m.group(1)
            mw_results[potential][limit]["estimate"] = estimate
            print("ESTIMATE: {}".format(estimate))
            conf_int = line4.split('[1] ')[1]
            mw_results[potential][limit]["conf_interval"] = conf_int
            print("CONFIDENCE INTERVAL: {}".format(conf_int))
            out.close()

            rm_temp = 'rm temp/temp.R'
            os.system(rm_temp)
            rm_out = 'rm temp/temp.out'
            os.system(rm_out)

    return mw_results


def calculate_mann_whitney_2_rows(filtered):
    '''
    Calculates the Mann Whitney test for the pairs of Native - All decoys values
    '''

    mw_results = {}

    for potential in filtered:

        mw_results.setdefault(potential, {})

        for limit in ('2_8', 'more'):

            mw_results[potential].setdefault(limit, {})

            native_list = filtered[potential]['native'][limit]['pred']
            docking_list = filtered[potential]['docking'][limit]['pred']
            temp_file = 'temp/temp.R'
            temp = open(temp_file, 'w')

            for values_list, type_list in ((native_list, 'native'), (docking_list, 'docking')):
                temp.write('{} <- c('.format(type_list))
                for x in xrange(len(values_list)):
                    if x == len(values_list) - 1:
                        temp.write('{})\n'.format(values_list[x]))
                    else:
                        temp.write('{}, '.format(values_list[x]))
            temp.write('mwt <- wilcox.test(native, docking, conf.int = TRUE)\n')
            temp.write('mwt$p.value\n')
            temp.write('mwt$estimate\n')
            temp.write('mwt$conf.int[1:2]\n')
            temp.close()

            command = '/soft/R/R-3.3.0/bin/Rscript temp/temp.R > temp/temp.out'
            os.system(command)

            out = open('temp/temp.out', 'r')
            num_line = 1
            for line in out:
                if num_line == 1:
                    line1 = line.strip()
                if num_line == 3:
                    line3 = line.strip()
                if num_line == 4:
                    line4 = line.strip()
                num_line+=1
            pval = line1.split('[1] ')[1]
            mw_results[potential][limit]["pvalue"] = pval
            print("PVAL: {}".format(pval))
            estimate = ""
            estimate_pattern = re.compile('\s*([-+]?[0-9]*\.?[0-9]*)')
            m = estimate_pattern.search(line3)
            if m:
                estimate = m.group(1)
            mw_results[potential][limit]["estimate"] = estimate
            print("ESTIMATE: {}".format(estimate))
            conf_int = line4.split('[1] ')[1]
            mw_results[potential][limit]["conf_interval"] = conf_int
            print("CONFIDENCE INTERVAL: {}".format(conf_int))
            out.close()

            rm_temp = 'rm temp/temp.R'
            os.system(rm_temp)
            rm_out = 'rm temp/temp.out'
            os.system(rm_out)

    return mw_results


def print_mann_whitney(mw_results_T, mw_results_R, mw_results_F, outfile):
    '''
    Prints the table of the Mann Whitney test results
    '''

    for potential in ['FiberDock', 'aVdW', 'rVdW', 'aElec', 'rElec', 'laElec', 'lrElec', 'HB', 'EPAIR', 'ES3DC', 'E3D']:

        potential_name = potential
        if potential == 'ROSSETADOCK':
            potential_name = 'RosettaDock'

        for percent in ['1.4', '2.8', '4.2', '>4.2']:
            if percent == '1.4':
                outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format( potential_name, '<=1.4', mw_results_T[potential]['1_4']["pvalue"], mw_results_T[potential]['1_4']["estimate"], mw_results_T[potential]['1_4']["conf_interval"],
                                                                                                                        mw_results_R[potential]['1_4']["pvalue"], mw_results_R[potential]['1_4']["estimate"], mw_results_R[potential]['1_4']["conf_interval"],
                                                                                                                        mw_results_F[potential]['1_4']["pvalue"], mw_results_F[potential]['1_4']["estimate"], mw_results_F[potential]['1_4']["conf_interval"] ))
            elif percent == '2.8':
                outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format( '', '<=2.8', mw_results_T[potential]['2_8']["pvalue"], mw_results_T[potential]['2_8']["estimate"], mw_results_T[potential]['2_8']["conf_interval"],
                                                                                                                        mw_results_R[potential]['2_8']["pvalue"], mw_results_R[potential]['2_8']["estimate"], mw_results_R[potential]['2_8']["conf_interval"],
                                                                                                                        mw_results_F[potential]['2_8']["pvalue"], mw_results_F[potential]['2_8']["estimate"], mw_results_F[potential]['2_8']["conf_interval"] ))
            elif percent == '4.2':
                outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format( '', '<=4.2', mw_results_T[potential]['4_2']["pvalue"], mw_results_T[potential]['4_2']["estimate"], mw_results_T[potential]['4_2']["conf_interval"],
                                                                                                                        mw_results_R[potential]['4_2']["pvalue"], mw_results_R[potential]['4_2']["estimate"], mw_results_R[potential]['4_2']["conf_interval"],
                                                                                                                        mw_results_F[potential]['4_2']["pvalue"], mw_results_F[potential]['4_2']["estimate"], mw_results_F[potential]['4_2']["conf_interval"] ))
            elif percent == '>4.2':
                outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format( '', '>4.2', mw_results_T[potential]['more']["pvalue"], mw_results_T[potential]['more']["estimate"], mw_results_T[potential]['more']["conf_interval"],
                                                                                                                        mw_results_R[potential]['more']["pvalue"], mw_results_R[potential]['more']["estimate"], mw_results_R[potential]['more']["conf_interval"],
                                                                                                                        mw_results_F[potential]['more']["pvalue"], mw_results_F[potential]['more']["estimate"], mw_results_F[potential]['more']["conf_interval"] ))

    return


def print_mann_whitney_2_rows(mw_results_T, mw_results_R, mw_results_F, outfile):
    '''
    Prints the table of the Mann Whitney test results
    '''

    for potential in ['FiberDock', 'aVdW', 'rVdW', 'aElec', 'rElec', 'laElec', 'lrElec', 'HB', 'EPAIR', 'ES3DC', 'E3D']:

        potential_name = potential
        if potential == 'ROSSETADOCK':
            potential_name = 'RosettaDock'

        for percent in ['2.8', '>2.8']:
            if percent == '2.8':
                outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format( potential_name, '<=2.8', mw_results_T[potential]['2_8']["pvalue"], mw_results_T[potential]['2_8']["estimate"], mw_results_T[potential]['2_8']["conf_interval"],
                                                                                                                        mw_results_R[potential]['2_8']["pvalue"], mw_results_R[potential]['2_8']["estimate"], mw_results_R[potential]['2_8']["conf_interval"],
                                                                                                                        mw_results_F[potential]['2_8']["pvalue"], mw_results_F[potential]['2_8']["estimate"], mw_results_F[potential]['2_8']["conf_interval"] ))
            elif percent == '>2.8':
                outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format( potential_name, '>2.8', mw_results_T[potential]['more']["pvalue"], mw_results_T[potential]['more']["estimate"], mw_results_T[potential]['more']["conf_interval"],
                                                                                                                        mw_results_R[potential]['more']["pvalue"], mw_results_R[potential]['more']["estimate"], mw_results_R[potential]['more']["conf_interval"],
                                                                                                                        mw_results_F[potential]['more']["pvalue"], mw_results_F[potential]['more']["estimate"], mw_results_F[potential]['more']["conf_interval"] ))

    return


if __name__ == "__main__":
    main()