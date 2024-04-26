# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
January 5th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

Code inspired by:
    https://nanohub.org/resources/struc2lammpsdf
"""


##############################
# Import Necessary Libraries #
##############################
import src.atom_typing.Gasteiger.read_parameters as read_parameters
import os

# Function to compute Gasteiger charges
def compute_charges(mm, chargefile, log):
    # Print statement
    log.out('\n\n\n')
    log.out('Starting Gasteiger algorithm to find charges. WARNING the algortihm itself is independent of atomic')
    log.out('positions and/or geometries, but the look up tables where Gasteirger parameters are intialized')
    log.out('from requires knowledge of hybridization of certain elements (C, N, O, ...). The method employeed')
    log.out('in this code is a hybrid approach by using bonding-topology (indepenent of atomic positions) and')
    log.out('a minimization of VSEPR ideal angles (dependent on atomic positions) onto local angles of each')
    log.out('atom in the system to find the hybridized state. This means for the best results of computed charge')
    log.out('that the geometries (set by bonding and atomic positions) of your file should be consistant with the')
    log.out('hybridized state of each atom: (IE Sp1 avg-angle ~180, Sp2 avg-angle ~120, Sp3 avg-angle ~109.5).')
    log.out('When bonding-topology could not determine the hybridized state the VSEPR ideal angle difference')
    log.out('from the avg-angle of bonded neighbors in the system will be minimized, such that there will always')
    log.out('be a predicted hybridized state.')
    log.out('')
    
    
    # Read parameters from file
    if os.path.isfile(chargefile):
        parameters = read_parameters.read(chargefile)
        log.out(f'Read in Gasteiger parameters from {chargefile}')
    else: log.error(f'ERROR Gasteiger parameter file: {chargefile} does not exist')
    
    
    # Initialize (sorted contigous atomIDs and lst)
    charge = {}; qn = {}; x= {}; xp = {}; # { atomID : value}
    sorted_contiguous_atomIDs = sorted([i for i in mm.atoms])
    IDs2skip = {i:False for i in sorted_contiguous_atomIDs}
    for ID in sorted_contiguous_atomIDs:
        
        # Set intialize as zeros and udpate if found
        ai=0; bi=0; ci=0;
        
        # Find atom element and hybridzation state
        hybridization = mm.atoms[ID].hybridization
        element = mm.atoms[ID].element
        
        # Try getting parameters; except warn user
        try:
            parms = parameters[element][hybridization]
            ai = parms.a; bi = parms.b; ci = parms.c;
        except:
            IDs2skip[ID] = True
            log.warn(f'WARNING atomID: {ID}   element: {element}   hybridization: {hybridization}   Gasteiger parameters could not be found and set to zeros')
        
        # Log values
        charge[ID] = 0; qn[ID] = 0; x[ID] = ai; xp[ID] = [ai, bi, ci];
        
        
    # Start iterative alogorithm
    max_iter = 500 # set max iteration
    log.out('\nFinding Gasteiger Charge ....')
    for n in range(1, max_iter):
        for ID in mm.bonds:
            id1, id2 = mm.bonds[ID].atomids
            id1_element = mm.atoms[id1].element
            id2_element = mm.atoms[id2].element
            
            d1=xp[id1][0] + xp[id1][1] + xp[id1][2]
            d2=xp[id2][0] + xp[id2][1] + xp[id2][2]
            d=min(d1, d2)
            if id1_element == 'H' or id2_element == 'H': d=max(d1, d2)
            if d==0: d=1e-6
            
            # Only change qn if atomID has Gasteiger parameters
            if not IDs2skip[id1] and not IDs2skip[id2]:
                qn[id1]=qn[id1] + ((x[id2] - x[id1])/d)*(0.5)**n
                qn[id2]=qn[id2] + ((x[id1] - x[id2])/d)*(0.5)**n

            
        converged = 0
        for i in sorted_contiguous_atomIDs:
            if abs(qn[i]) < 1e-8: converged += 1
            charge[i] = charge[i] + qn[i]
            qn[i] = 0.0
            x[i]=xp[i][0] + xp[i][1]*charge[i] + xp[i][2]*charge[i]*charge[i]
        
        # Break out condition
        if converged == len(sorted_contiguous_atomIDs):
            log.out('\nGasteiger method converged'); break
        
        
    # Check for net-zero charge; if not scale charges
    total_charge = sum(list(charge.values()))
    
    # IF total charge is not close to zero, scale all charges
    if abs(total_charge) > 1e-10:
        log.out('Total charge of the system before charge rescale: {:^4.6f}'.format( round(total_charge, 6)))
        for i in charge:
            charge[i] -= total_charge/len(charge)
        log.out('Total charge of the system after charge rescale: {:^4.6f}'.format(round(sum(list(charge.values())), 6)))
    else:
        log.out('Total charge of the system: {:^4.6f}'.format(round(total_charge,6)))
        
        
    # Print Citations to users  
    log.out('')
    log.out('Cite charge method:')
    log.out('    Ref : J.Gasteiger, M. Marseli, Iterative Equalization of Oribital Electronegatiity A Rapid')
    log.out('          Access to Atomic Charges, Tetrahedron Vol 36 p3219 1980')
    
    # Update mm.atoms[atomID].charge
    for i in mm.atoms:
        mm.atoms[i].charge = charge[i]            
    return mm