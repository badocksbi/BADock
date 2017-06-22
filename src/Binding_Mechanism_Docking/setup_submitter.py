from BioLib.Tools.Submitters import *
import optparse, os

def parse_options():
    '''
    Parse arguments
    @return options object
    '''
    parser = optparse.OptionParser('setup_submitter.py -s SCRIPT_NAME [-n] [-g]')
    parser.add_option("-s", dest="script_name", action="store", choices=['decoys', 'location', 'score'], help="The script to submit [MANDATORY]", metavar="SCRIPT_NAME")
    parser.add_option("-n", dest="nips", action="store_true", default=False, help="Non-interacting proteins setup")
    parser.add_option("-g", dest="gaudi", action="store_true", default=False, help="Send Jobs to gaudi (must be in node master)")

    (options, args) = parser.parse_args()
    return options

#--------------------------#
#          MAIN            #
#--------------------------#

pdb_path = os.path.expanduser('~/Databases/DockingBenchmark5')
ppi_file = os.path.join(pdb_path, 'DBM5.txt')
jobs_path = os.path.abspath('benchmark5_jobs')
log_path = os.path.abspath('logs')

options = parse_options()
if options.nips:
    pdb_path += '_randoms'
    jobs_path += '_nips'
    log_path += '_nips'
if options.gaudi:
    submitter = GaudiSubmitter(qsub='bigmem', log_path=log_path)
else:
    submitter = LocalSubmitter(log_path=log_path)

ppi_fo = open(ppi_file, 'r')
for ppi in ppi_fo:
    ppi = ppi.strip('\n').split('\t')[0][0:4]
    job_path = os.path.join(jobs_path, ppi)
    if not os.path.isdir(job_path):
        os.mkdir(job_path)
        os.mkdir(job_path+'/'+'decoys')
    if options.script_name == 'decoys':
        cmd = '/soft/devel/python-2.7/bin/python %s -i %s -j %s -p %s' % (os.path.join(os.path.dirname(__file__), 'setup_decoys.py'), ppi, job_path, pdb_path)
        job_id = 'job.%s.setup_decoys' % ppi
    if options.script_name == 'location':
        cmd = '/soft/devel/python-2.7/bin/python %s -i %s -j %s -p %s' % (os.path.join(os.path.dirname(__file__), 'setup_location.py'), ppi, job_path, pdb_path)
        job_id = 'job.%s.setup_location' % ppi
    if options.script_name == 'iLoops':
        cmd = '/soft/devel/python-2.7/bin/python %s -i %s -j %s -p %s' % (os.path.join(os.path.dirname(__file__), 'setup_score.py'), ppi, job_path, pdb_path)
        job_id = 'job.%s.setup_iLoops' % ppi
    submitter.submit(cmd, job_id)
ppi_fo.close()
