from BioLib.Tools.Submitters import *
import optparse, os

def parse_options():
    '''
    Parse arguments
    @return options object
    '''
    parser = optparse.OptionParser('setup_submitter.py -s SCRIPT_NAME [-g]')
    parser.add_option("-s", dest="script_name", action="store", choices=['docking', 'score'], help="The script to submit [MANDATORY]", metavar="SCRIPT_NAME")
    parser.add_option("-g", dest="gaudi", action="store_true", default=False, help="Send Jobs to gaudi (must be in node master)")
    (options, args) = parser.parse_args()

    if options.script_name == None:
        parser.error('Missing arguments\n')
    return options

#--------------------------#
#          MAIN            #
#--------------------------#

pdbs_path = os.path.expanduser('~/Databases/AffinityBenchmark')
ba_file = os.path.join(pdbs_path, 'affinity_benchmark_all_dG.tsv')
jobs_path = os.path.abspath('ba_jobs')
log_path = os.path.abspath('logs')

options = parse_options()
if options.gaudi:
    submitter = GaudiSubmitter(qsub='sbi', log_path=log_path)
else:
    submitter = LocalSubmitter(log_path=log_path)

bonvin_cleaned_dataset = get_bonvin_cleaned_dataset()

ba_fo = open(ba_file, 'r')
for ppi in ba_fo:
    ppi = ppi.strip('\n').split('\t')[0]
    job_path = os.path.join(jobs_path, ppi)
    if not os.path.isdir(job_path):
        os.mkdir(job_path)
    if options.script_name == 'docking':
        cmd = '/soft/devel/python-2.7/bin/python %s -i %s -j %s -p %s' % (os.path.join(os.path.dirname(__file__), 'setup_docking.py'), ppi, job_path, pdbs_path)
        job_id = 'job.%s.setup_docking' % ppi
    if options.script_name == 'score':
        cmd = '/soft/devel/python-2.7/bin/python %s -i %s -j %s -p %s' % (os.path.join(os.path.dirname(__file__), 'setup_score.py'), ppi, job_path, pdbs_path)
        job_id = 'job.%s.setup_score' % ppi
    submitter.submit(cmd, job_id)
ba_fo.close()
