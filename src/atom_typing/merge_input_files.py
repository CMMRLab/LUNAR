# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
January 5th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.bonds_via_distance as bonds_via_distance
import src.mol2SYBYL2lmp as mol2SYBYL2lmp
import src.read_reaxff as read_reaxff
import src.read_lmp as read_lmp
import src.mol2lmp as mol2lmp
import src.pdb2lmp as pdb2lmp
import os

# Try importing rdkit addins
try:
    import src.rdkit_optional_addins as rdkit_optional_addins
    rdkit_flag = True
except: rdkit_flag = False




# Function to read -> merge -> add missing data from each file that atom_typing will be able to read
def merge(topofile, bondfile, mass_map, bondorder, maxbonded, boundary, vdw_radius_scale, reset_charges, bonds_via_distance_override, log):
        
    ################################################
    # read topofile and bondfile (when applicable) #
    ################################################
    reaxff_flag = False # will update if files are for reaxff conversion
    bonddist_flag = False # will update if files have no bonds and bonds are found via distance
    # Read lammps data file
    if topofile.endswith('data') or topofile.endswith('dat') or  topofile.endswith('data.gz') or topofile.endswith('dat.gz'):
        if os.path.isfile(topofile):
            m = read_lmp.Molecule_File(topofile, method='forward', sections = ['Atoms', 'Bonds'])
            log.out(f'Read in {m.filename} LAMMPS datafile')
        else: log.error(f'ERROR lammps datafile: {topofile} does not exist')
            
        # If read in a LAMMPS .data file it will not have an element 
        # attribute in m.atoms so find element symbol from mass_map
        log.out('\n\nUsing mass_map dictionary to set element symbols ...')
        for i in m.atoms:
            # Find element from mass_map and add to atom instance
            atom = m.atoms[i]; mass = m.masses[atom.type].coeffs[0];
            try:
                atom.element = [i for i in mass_map if mass in mass_map[i]][0];
                atom.comment = [i for i in mass_map if mass in mass_map[i]][0];
            except: log.error(f'ERROR Not all masses in {topofile} are in the mass_map dictionary. Failed for mass: {mass}')
            
        # Add element symbol to m.masses[ID].type
        for i in m.masses:
            mass = m.masses[i]
            mass.type = [i for i in mass_map if mass.coeffs[0] in mass_map[i]][0]
            
        # Check if len(m.bonds) == 0, if so it is a reaxff file and read
        # the reaxff orginal bond order file and add to m.bonds instance
        if len(m.bonds) == 0 and bondfile != 'n.u.':
            reaxff_flag = True
            
            # update bondfile name if bondfile starts with topofile build bondfile actual name
            if bondfile.startswith('topofile'):            
                base = topofile[:topofile.rfind('.')]
                ext = bondfile.split('.')[-1]
                bondfile = '{}.{}'.format(base, ext)
                log.out(f'Using path and topfile name from topofile to set bondfile = {bondfile}')
 
            
            # Generate bond_info dictionary
            bond_info = {} # {atomtype: ('element symbol', max number of bonded atoms)} 
            for i in m.masses:
                mass = m.masses[i].coeffs[0]; element = [i for i in mass_map if mass in mass_map[i]][0];
                bond_info[i] = (element, maxbonded[element])

            # Find bonding information from reaxff file and abos/nlps avgs class
            reaxff = read_reaxff.create_bonds(bondfile, bond_info, bondorder, log, maxstep='all')
            log.out(f'Read in {bondfile} LAMMPS reaxff bond-order file')

            
            # Add bonds to m.bonds and update m.nbonds and m.nbondtypes
            class Bonds: pass # .type .atomids
            for n, i in enumerate(reaxff.bonds):
                b = Bonds()
                b.type = 1 # set as 1 for now ...
                b.atomids = sorted(i)
                m.bonds[n+1] = b
            m.nbonds = len(m.bonds); m.nbondtypes = 1;
            
            # Add avgs, abo_stats, statistics to reaxff class and append to m
            class reaxff_stats(): pass
            R = reaxff_stats()
            R.abo_stats = reaxff.abo_stats
            R.bo_stats = reaxff.statistics
            R.ntimesteps = len(reaxff.timesteps)
            R.nbonds = len(reaxff.bonds)
            R.nflaggedbonds = len(reaxff.flagged_bonds)
            R.maxbonded = maxbonded
            m.reaxff = R

    # Read .mol file
    elif topofile.endswith('mol') or topofile.endswith('sdf'):
        if os.path.isfile(topofile):
            m = mol2lmp.Molecule_File(topofile)
            log.out(f'Read in {m.filename} chemdraw .mol or .sdf file')
        else: log.error(f'ERROR .mol or .sdf file: {topofile} does not exist')
            
    # Read .mol2 file (VMD MDL MOL2 file)
    elif topofile.endswith('mol2'):
        if os.path.isfile(topofile):
            m = mol2SYBYL2lmp.Molecule_File(topofile)
            log.out(f'Read in {m.filename} SYBYL MOL2 file')
        else: log.error(f'ERROR .mol2 file: {topofile} does not exist')
            
    # Read .pdb file
    elif topofile.endswith('pdb'):
        if os.path.isfile(topofile):
            m = pdb2lmp.Molecule_File(topofile)
            log.out(f'Read in {m.filename} pdb file')
        else: log.error(f'ERROR .pdb file: {topofile} does not exist')
            
    # Read .smiles string only if rdkit was determined to be installed
    elif topofile.endswith('smiles'):
        if rdkit_flag:
            mol = rdkit_optional_addins.rdkit2lmp(topofile, addHs=True, gen_mclass=True)
            m = mol.m
            log.out(f'Read in {topofile} SMILES string. Cite Rdkit for functionality:')
            log.out('    RDKit: Open-source cheminformatics; http://www.rdkit.org')
        else: log.error('ERROR Requesting rdkit addin and rdkit is not currently installed')
            
            
    # If topofile was a .mol or .sdf or .mol2 file it will not have m.masses nor correct atom types so find elements and set types
    if topofile.endswith('mol') or topofile.endswith('sdf') or topofile.endswith('mol2') or topofile.endswith('pdb') or topofile.endswith('smiles'):
        # Find all elements in system to set atomtypes
        elements = sorted({m.atoms[i].element for i in m.atoms})
        atomtypes = {i:n+1 for n, i in enumerate(elements)}

        # Set mass from for each atomtype in atomtypes and
        # build m.masses instance as read_lmp would have
        class Coeff_class: pass  # .type .coeffs = []
        m.masses = {}
        for i in atomtypes:
            c = Coeff_class()
            c.type = i # set type as element
            c.coeffs = [mass_map[i][0]] # set mass mass_map
            m.masses[atomtypes[i]] = c
        m.natomtypes = len(m.masses)
        
        # Assign atom type in m.atoms based on element
        for i in m.atoms:
            atom = m.atoms[i];
            atom.type = atomtypes[atom.element];
            atom.comment = atom.element
            
    # Check if len(m.bonds) == 0, if so generate bonds via vdw distance cutoffs
    if len(m.bonds) == 0 or bonds_via_distance_override:
        log.out('Read in file does not contain bonds. Bonds will be found via distance searching using vdw radius.')
        log.out('vdw_radius_scale in inputs will set scale value of distances to search for bond determination.')
        log.out('maxbonded will then be used to reduce the number of bonds down to an acceptable amount if to many')
        log.out('atoms are bonded to a given element type. The bonded atoms that are furthest away will be rejected.')
        bonddist_flag = True
        
        # Finding bonds
        bond_creation = bonds_via_distance.generate(m, boundary, vdw_radius_scale, maxbonded, log)
        
        # Add bonds to m.bonds and update m.nbonds and m.nbondtypes
        class Bonds: pass # .type .atomids
        for n, i in enumerate(bond_creation.bonds):
            b = Bonds()
            b.type = 1 # set as 1 for now ...
            b.atomids = sorted(i)
            m.bonds[n+1] = b
        m.nbonds = len(m.bonds); m.nbondtypes = 1;
        
        # Add avgs, abo_stats, statistics to reaxff class and append to m
        class bonddist_stats(): pass
        b = bonddist_stats()
        b.dist_stat = bond_creation.statistics
        b.nbonds = len(bond_creation.bonds)
        b.nflaggedbonds = len(bond_creation.flagged_bonds)
        b.maxbonded = maxbonded
        b.boundary = boundary
        b.images = bond_creation.images
        b.nb_count = bond_creation.nb_count
        b.maxbond = bond_creation.maxbond
        b.bond_status = bond_creation.bond_status
        b.vdw_radius_scale = vdw_radius_scale
        m.bonds_via_dist = b
            
    # create a new instance in m called reaff_flag to keep track of typing inputs
    m.reaxff_flag = reaxff_flag
    
    # create a new instance in m called bonddist_flag to keep track of input creation
    m.bonddist_flag = bonddist_flag
    
    # find all elements in file and add to m
    m.elements = sorted({m.masses[i].type for i in m.masses})         
    return m
