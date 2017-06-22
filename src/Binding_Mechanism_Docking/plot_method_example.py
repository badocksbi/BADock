from BioLib import *
import optparse, os, sys

def get_rmsd(job_path, pdb_name):
    '''
    Get the RMSD for ach decoy using the rmsd_splitpotentials file
    '''
    # Check if needed files exists #
    rmsd_path = os.path.join(job_path, pdb_name, 'rmsd_splitpotentials.txt')
    if not os.path.isfile(rmsd_path):
        raise RuntimeError('NO RMSD FILE FOUND')
        exit()
    # Get rmsd #
    rmsd_decoys = {}
    rmsd_decoys_fo = open(rmsd_path, 'r')
    for decoy in rmsd_decoys_fo:
        decoy = decoy.strip().split('\t')
        rmsd_decoys[decoy[0]] = float(decoy[1])
    rmsd_decoys_fo.close()
    return rmsd_decoys

def get_decoys_byclass(job_path, pdb_name):
    '''
    Get the decoys divided by class
    '''
    # Check if needed files exists #
    decoys_location_path = os.path.join(job_path, pdb_name, 'decoys_location.txt')
    if not os.path.isfile(decoys_location_path):
        raise RuntimeError('NO DECOYS LOCATION FILE FOUND')
        exit()
    # Get rmsds #
    location_decoys = {'N': [],'I': [], 'P': [], 'E': []}
    rmsd_decoys = get_rmsd(job_path, pdb_name)
    # Get location decoys #
    decoys_locations_fo = open(decoys_location_path, 'r')
    for decoy in decoys_locations_fo:
        decoy = decoy.strip('\n').split('\t')
        if not decoy[0] in rmsd_decoys:
            continue
        if rmsd_decoys[decoy[0]] < 10:
            location_decoys['N'].append(int(decoy[0]))
        else:
            location_decoys[decoy[1]].append(int(decoy[0]))
    decoys_locations_fo.close()
    return location_decoys

def get_structures(rec_path, lig_path):
    '''
    Parse the original PDB files in order to get the structure objects
    '''
    try:
        rec_struc = PDB.read_pdb(rec_path, merge_chains=True)
        rec_struc.clean()
        lig_struc = PDB.read_pdb(lig_path, merge_chains=True)
        lig_struc.clean()
    except Exception as e:
        raise RuntimeError('ERROR: Cannot parse PDB: %s\n' % str(e))
        exit()
    return rec_struc, lig_struc

def get_patchdock(job_path, pdb_name):
    '''
    PARSE PATCHDOCK transformation file
    '''
    patchdock_path = os.path.join(job_path, pdb_name, 'patchdock_transform.txt')
    if not os.path.isfile(patchdock_path):
        raise RuntimeError('ERROR: No PATCHDOCK TRASNFORMATION file found\n')
        exit()
    try:
        patchdock_results = PATCHDock(patchdock_path)
    except Exception as e:
        raise RuntimeError('ERROR: Cannot parse PATCHDOCK TRASNFORMATION file: %s\n' % str(e))
        exit()
    return patchdock_results

#--------------------------#
#          MAIN            #
#--------------------------#

if __name__ == "__main__":

    # Settings #
    pdb_name = '1JTG'
    job_path = os.path.abspath('benchmark5_jobs')
    pdbs_path = os.path.expanduser('~/Databases/DockingBenchmark5')
    # Get the decoy numbers by class #
    decoys_byclass = get_decoys_byclass(job_path, pdb_name) 
    if not decoys_byclass['N'] or not decoys_byclass['I'] or not decoys_byclass['P'] or not decoys_byclass['E']:
        raise RuntimeError('NO NATIVE, INTERFACE, PARTIAL OR EXTERNAL DECOYS')
    # Get structures #
    rec_path = os.path.join(pdbs_path, pdb_name+'_r_b.pdb')
    lig_path = os.path.join(pdbs_path, pdb_name+'_l_b.pdb')
    rec_struc, lig_struc = get_structures(rec_path, lig_path)
    # Print first decoy of each class #
    patchdock = get_patchdock(job_path, pdb_name)
    patchdock.get_decoy(decoys_byclass['N'][0]).print_structure(rec_struc, lig_struc, pdb_name+'_native')
    patchdock.get_decoy(decoys_byclass['I'][0]).print_structure(rec_struc, lig_struc, pdb_name+'_interface')
    patchdock.get_decoy(decoys_byclass['P'][1]).print_structure(rec_struc, lig_struc, pdb_name+'_partial')
    patchdock.get_decoy(decoys_byclass['E'][0]).print_structure(rec_struc, lig_struc, pdb_name+'_external')
    # Print interacting residues #
    ppi_struc = Interaction(rec_struc, lig_struc)
    int_residues = ppi_struc.get_interacting_residues()
    sys.stdout.write('color purple :')
    for i in int_residues[0]:
        sys.stdout.write('%d.%s,' % (i.get_num(), rec_struc.get_chain()))
    sys.stdout.write('\ncolor cornflower blue :')
    for i in int_residues[1]:
        sys.stdout.write('%d.%s,' % (i.get_num(), lig_struc.get_chain()))
    sys.stdout.write('\n')





