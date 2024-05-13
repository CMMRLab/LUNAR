# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.3
April 1st, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
# Class to rebuild m class with atomIDs removed and remapped with atoms, bonds, angles, dihedrals, impropers, and coeffs
class Atom: pass # .type .molid .charge .x .y .z .ix. iy .iz .comment
class Bond: pass  # .type .atomids = [atom1id, atom2id]
class Angle: pass  # .type .atomids = [atom1id, atom2id, atom3id]
class Dihedral:pass  # .type .atomids = [atom1id, atom2id, atom3id, atom4id]
class Improper: pass  # .type .atomids = [atom1,atom2,atom3,atom4]
class Coeff_class: pass  # .type .coeffs = []
class constructor:
    def __init__(self, m, atoms2remove, log, method='atomIDs', pflag=True):
        # Structure info
        self.atoms = {}  # {atom number : atom object}
        self.bonds = {}  # {bond number : bond object}
        self.angles = {}  # {angle number : angle object}
        self.dihedrals = {}  # {dihedral number : dihedral object}
        self.impropers = {}  # {improper number : improper object}
        self.filename = m.filename
        
        # Quantities
        self.natoms = 0
        self.nbonds = 0
        self.nangles = 0
        self.ndihedrals = 0
        self.nimpropers = 0
        self.natomtypes = m.natomtypes
        self.nbondtypes = m.nbondtypes
        self.nangletypes = m.nangletypes
        self.ndihedraltypes = m.ndihedraltypes
        self.nimpropertypes = m.nimpropertypes
        self.nbondbond = m.nbondbond
        self.nbondangle = m.nbondangle
        self.nangleangle = m.nangleangle
        self.nangleangletorsion = m.nangleangletorsion
        self.nendbondtorsion = m.nendbondtorsion
        self.nmiddlebondtorsion = m.nmiddlebondtorsion
        self.nbondbond13 = m.nbondbond13
        self.nangletorsion = m.nangletorsion
        self.xbox_line = m.xbox_line
        self.ybox_line = m.ybox_line
        self.zbox_line = m.zbox_line
        self.xy = m.xy
        self.xz = m.xz
        self.yz = m.yz
        
        # Coeffs
        self.header = m.header
        self.masses = m.masses; self.mass_coeffs_style_hint = m.mass_coeffs_style_hint;
        self.pair_coeffs = m.pair_coeffs; self.pair_coeffs_style_hint = m.pair_coeffs_style_hint;
        self.bond_coeffs = m.bond_coeffs; self.bond_coeffs_style_hint = m.bond_coeffs_style_hint;
        self.angle_coeffs = m.angle_coeffs; self.angle_coeffs_style_hint = m.angle_coeffs_style_hint;
        self.dihedral_coeffs = m.dihedral_coeffs; self.dihedral_coeffs_style_hint = m.dihedral_coeffs_style_hint;
        self.improper_coeffs = m.improper_coeffs; self.improper_coeffs_style_hint = m.improper_coeffs_style_hint;
        self.bondbond_coeffs = m.bondbond_coeffs; self.bondbond_coeffs_style_hint = m.bondbond_coeffs_style_hint;
        self.bondangle_coeffs = m.bondangle_coeffs; self.bondangle_coeffs_style_hint = m.bondangle_coeffs_style_hint;
        self.angleangletorsion_coeffs = m.angleangletorsion_coeffs; self.angleangletorsion_coeffs_style_hint = m.angleangletorsion_coeffs_style_hint;
        self.endbondtorsion_coeffs = m.endbondtorsion_coeffs; self.endbondtorsion_coeffs_style_hint = m.endbondtorsion_coeffs_style_hint;
        self.middlebondtorsion_coeffs = m.middlebondtorsion_coeffs; self.middlebondtorsion_coeffs_style_hint = m.middlebondtorsion_coeffs_style_hint;
        self.bondbond13_coeffs = m.bondbond13_coeffs; self.bondbond13_coeffs_style_hint = m.bondbond13_coeffs_style_hint;
        self.angletorsion_coeffs = m.angletorsion_coeffs; self.angletorsion_coeffs_style_hint = m.angletorsion_coeffs_style_hint;
        self.angleangle_coeffs = m.angleangle_coeffs; self.angleangle_coeffs_style_hint = m.angleangle_coeffs_style_hint;
        self.atom_type_labels_forward = m.atom_type_labels_forward; self.atom_type_labels_reverse = m.atom_type_labels_reverse;
        self.bond_type_labels_forward = m.bond_type_labels_forward; self.bond_type_labels_reverse = m.bond_type_labels_reverse;
        self.angle_type_labels_forward = m.angle_type_labels_forward; self.angle_type_labels_reverse = m.angle_type_labels_reverse;
        self.dihedral_type_labels_forward = m.dihedral_type_labels_forward; self.dihedral_type_labels_reverse = m.dihedral_type_labels_reverse;
        self.improper_type_labels_forward = m.improper_type_labels_forward; self.improper_type_labels_reverse = m.improper_type_labels_reverse;
        
        # Print method of identifying atoms
        if pflag:
            if method == 'atomIDs':
                log.out('    Identifying atoms based on atomIDs')
            elif method == 'typeIDs':
                log.out('    Identifying atoms based on TypeIDs')
            elif method == 'cluster-mass':
                log.out('    Identifying atoms based on mass of cluster')
            elif method == 'cluster-size':
                log.out('    Identifying atoms based on size of cluster')
            else:
                log.error(f"ERROR trying to remove atoms with method = {method}, currently supported methods are 'atomIDs' or 'typeIDs'")
                
        # If method is mass perform cluster analsis
        if method in ['cluster-mass', 'cluster-size']:
            # Check for proper inputs (only one value)
            if len(atoms2remove) != 1:
                log.error(f"ERROR using {method} method, but input in atoms2remove {str(atoms2remove)} is not a single value.")
            elif len(m.bonds) == 0:
                log.error(f"ERROR using {method} method, no bonds exists to perform cluster analysis.")
            else:
                # Generate graph
                graph = {i:[] for i in m.atoms}
                for i in m.bonds:
                    id1, id2 = m.bonds[i].atomids
                    graph[id1].append(id2)
                    graph[id2].append(id1)
                
                # Find clusters
                clusters = set([]) # { (tuple of atoms in cluster1), (nclusters) }
                checked = {ID:False for ID in m.atoms}
                size_total = len(m.atoms); mass_total = 0
                for ID in graph:
                    mass_total += m.masses[m.atoms[ID].type].coeffs[0]
                    if checked[ID]: continue
                    visited=set([ID]); queue=[ID];
                    while queue:
                        s = queue.pop(0) 
                        for neighbor in graph[s]:
                            if checked[neighbor]: continue
                            visited.add(neighbor)
                            queue.append(neighbor)
                            checked[neighbor]=True
                    clusters.add( tuple(sorted(visited)) )
                clusters = sorted(clusters, key=lambda x: x[0]) # Sort all clusters based on 1st atomID in cluster
                clusters = sorted(clusters, key=len, reverse=True) # Sort all clusters by number of atoms
                
                # analyze clusters and add cluster_mass and cluster_size to atoms object
                def getmass(c):                                                                                                                  
                    return sum([m.masses[m.atoms[i].type].coeffs[0] for i in c])
                log.out('\n\n    ------------------------------------Cluster Analysis------------------------------------')
                log.out('    {:^10} {:^15} {:^15} {:^15} {:^15} {:^15}'.format('molID', 'size', '%size', 'mass', '%mass', 'status'))
                log.out('    ----------------------------------------------------------------------------------------')  
                for n, cluster in enumerate(clusters, 1):
                    # print table
                    psize = 0; pmass = 0; # Intialize as zeros and update
                    mass = getmass(cluster); size = len(cluster);
                    if size_total > 0: psize = 100*size/size_total
                    if mass_total > 0: pmass = 100*mass/mass_total
                    status = 'kept'
                    if method == 'cluster-mass' and mass <= atoms2remove[0]: status = 'removed'
                    if method == 'cluster-size' and size <= atoms2remove[0]: status = 'removed'
                    log.out('    {:^10} {:^15} {:^15.2f} {:^15.2f} {:^15.2f} {:^15}'.format(n, size, psize, mass, pmass, status))
                    
                    # add to m.atoms
                    for i in cluster:
                        m.atoms[i].cluster_mass = mass
                        m.atoms[i].cluster_size = size
                        m.atoms[i].cluster_ID = n
                log.out('\n\n')
            
        # Start removing atoms
        atomID_map = {} # { orginal atomID : new atomID }
        atoms = sorted(m.atoms.keys()) # sort atoms to keep atomIDs as close as possible
        flags = {i:False for i in atoms} # Flags to reduce bonds, angles, dihedrals, and impropers
        for i in atoms:
            saveatom = True
            if method == 'atomIDs':
                if i in atoms2remove: saveatom = False
            elif method == 'typeIDs':
                if m.atoms[i].type in atoms2remove: saveatom = False
            elif method == 'cluster-mass':
                if m.atoms[i].cluster_mass <= atoms2remove[0]: saveatom = False
            elif method == 'cluster-size':
                if m.atoms[i].cluster_size <= atoms2remove[0]: saveatom = False
            
            if saveatom:
                atom = m.atoms[i]
                self.natoms += 1 # increment atom count
                atomID_map[i] = self.natoms
                self.atoms[self.natoms] = atom
                flags[i] = True
        self.atomID_map = atomID_map # Creat atomID_map attribute
        if pflag: log.out('    Original number of atoms: {}    new number of atoms: {}'.format(len(m.atoms), len(self.atoms)))
        
        # Start removing bonds
        for i in m.bonds:
            bond = m.bonds[i]
            id1, id2 = bond.atomids
            if flags[id1] and flags[id2]:
                self.nbonds += 1
                b = Bond()
                b.type = bond.type
                b.atomids = [atomID_map[id1], atomID_map[id2]]
                self.bonds[self.nbonds] = b
        if pflag: log.out('    Original number of bonds: {}    new number of bonds: {}'.format(len(m.bonds), len(self.bonds)))
        
        # Start removing angles
        for i in m.angles:
            angle = m.angles[i]
            id1, id2, id3 = angle.atomids
            if flags[id1] and flags[id2] and flags[id3]:
                self.nangles += 1
                a = Angle()
                a.type = angle.type
                a.atomids = [atomID_map[id1], atomID_map[id2], atomID_map[id3]]
                self.angles[self.nangles] = a
        if pflag: log.out('    Original number of angles: {}    new number of angles: {}'.format(len(m.angles), len(self.angles)))
        
        # Start removing dihedrals
        for i in m.dihedrals:
            dihedral = m.dihedrals[i]
            id1, id2, id3, id4 = dihedral.atomids
            if flags[id1] and flags[id2] and flags[id3] and flags[id4]:
                self.ndihedrals += 1
                d = Dihedral()
                d.type = dihedral.type
                d.atomids = [atomID_map[id1], atomID_map[id2], atomID_map[id3], atomID_map[id4]]
                self.dihedrals[self.ndihedrals] = d
        if pflag: log.out('    Original number of dihedrals: {}    new number of dihedrals: {}'.format(len(m.dihedrals), len(self.dihedrals)))
        
        # Start removing impropers
        for i in m.impropers:
            improper = m.impropers[i]
            id1, id2, id3, id4 = improper.atomids
            if flags[id1] and flags[id2] and flags[id3] and flags[id4]:
                self.nimpropers += 1
                i = Improper()
                i.type = improper.type
                i.atomids = [atomID_map[id1], atomID_map[id2], atomID_map[id3], atomID_map[id4]]
                self.impropers[self.nimpropers] = i
        if pflag: log.out('    Original number of impropers: {}    new number of impropers: {}'.format(len(m.impropers), len(self.impropers)))