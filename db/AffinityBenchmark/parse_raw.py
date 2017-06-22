import os

ab_raw_fo = open('affinity_benchmark1_raw.tsv')
ab_dG_fo = open('affinity_benchmark1_dG2.tsv', 'w')
ab_raw_fo.readline()
for complex_ba in ab_raw_fo:
    complex_ba = complex_ba.strip('\n').split('\t')
    ab_dG_fo.write('%s\t%s\t%s\t%s\t%s\n' % (complex_ba[0][1:5], complex_ba[1][1:-1], complex_ba[7], complex_ba[8], complex_ba[9]))
ab_dG_fo.close()
ab_raw_fo.close()
