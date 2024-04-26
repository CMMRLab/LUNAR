# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 16th, 2022
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Function to write log file
def out(mm, bp, log, version, ff_name):
    ##########################################
    # Write the elements found in the system #
    ##########################################
    log.out('\n\nElements found in system:')
    for i in mm.elements:
        log.out('{} {}'.format('-', i))

    ####################################
    # If reaxff flag write reaxff info #
    ####################################
    if mm.reaxff_flag:
        # Write reaxff header
        log.out('\n\n\nReaxFF atom-typing specific section:')
        
        # Write intial reaxff parameters of nbonds, nflaggedbonds, and ntimesteps
        log.out('{} {} {}'.format('Each bonding atom ids have had', mm.reaxff.ntimesteps, 'Bond orders (BO) averaged over to find average BO for each bonding pair'))
        log.out('{} {} {}'.format('Total bonds created     :', mm.reaxff.nbonds, '(due to meeting specified criteria)'))
        log.out('{} {} {}\n'.format('Total bonds not created :', mm.reaxff.nflaggedbonds, '(due to not meeting specified criteria)'))
        
        # Write bond type stats table
        log.out('---------------------------Bond type bond order statistics and info---------------------------')
        log.out('{:^10} {:^10} {:^10} {:^15} {:^15} {:^15} {:^15}'.format('Bond', 'Bond', 'Bond Order', 'Bond Order', 'Bond Order', 'Bond Order','Cut-off'))
        log.out('{:^10} {:^10} {:^10} {:^15} {:^15} {:^15} {:^10}'.format('Type', 'Count', 'Average', 'Minimum', 'Maximum', 'Standard Deviation','used'))
        log.out('----------------------------------------------------------------------------------------------')
        for bond in mm.reaxff.bo_stats:
            i = mm.reaxff.bo_stats[bond]
            if i:
                log.out('{:^10} {:^10} {:^10} {:^15} {:^15} {:^15} {:^15}'.format(bond, i.count, i.avg, i.min, i.max, i.std, i.cutoff))
                 
        # Write abo type stats
        log.out('\n')
        log.out('-------------------------------Element abo statistics and info--------------------------------')
        log.out('{:<10} {:^8} {:^10} {:^15} {:^13} {:^15} {:^19}'.format('Element', 'Element', 'abo', 'abo', 'abo', 'abo','Max bonded'))
        log.out('{:<10} {:^8} {:^10} {:^15} {:^13} {:^15} {:^14}'.format('Type', 'Count', 'Average', 'Minimum', 'Maximum', 'Standard Deviation','Cut-off'))
        log.out('----------------------------------------------------------------------------------------------')
        for element in mm.reaxff.abo_stats:
            i = mm.reaxff.abo_stats[element]
            log.out('{:<10} {:^8} {:^10} {:^15} {:^13} {:^15} {:^19}'.format(element, i.count, i.avg, i.min, i.max, i.std, i.cutoff))
     
        
    ###################################################
    # If bond distance flag write bond distance stats #
    ###################################################
    if mm.bonddist_flag:
        # Write bond status and boundary used to create bonds
        log.out('\n{} {}'.format('vdw radius scale', mm.bonds_via_dist.vdw_radius_scale))
        log.out('{} {} {} {}'.format('boundary used:', mm.bonds_via_dist.boundary, '  nimages searched:', len(mm.bonds_via_dist.images)))
        log.out('{} {} {} {}'.format('non-periodic bonds found:', mm.bonds_via_dist.bond_status['non-periodic'], '    periodic bonds found:', mm.bonds_via_dist.bond_status['periodic']))

        
        # Write bond type stats table
        log.out('------------------------------------Bond type bond length statistics and info------------------------------------')
        log.out('{:^10} {:^10} {:^15} {:^20} {:^20} {:^20} {:^18}'.format('Bond', 'Bond', 'Bond Length', 'Bond Length', 'Bond Length', 'Bond Length','Cut-off'))
        log.out('{:^10} {:^10} {:^15} {:^20} {:^20} {:^20} {:^18}'.format('Type', 'Count', 'Average', 'Minimum', 'Maximum', 'Standard Deviation','used'))
        log.out('-----------------------------------------------------------------------------------------------------------------')
        for bond in mm.bonds_via_dist.dist_stat:
            i = mm.bonds_via_dist.dist_stat[bond]
            if i:
                log.out('{:^10} {:^10} {:^15} {:^20} {:^20} {:^20} {:^18}'.format(bond, i.count, i.avg, i.min, i.max, i.std, i.cutoff))
                
        # Write nb table
        columns = 'element'
        for i in range(mm.bonds_via_dist.maxbond+1):
            tmp = '{:>10} {:^5}'.format('nb', 'count')
            columns += '{:^16}: {:^3}'.format(tmp, i)
        log.out('\n\n{}'.format( ''.join(len(columns)*['-']) ))
        log.out('{}'.format(columns))
        log.out('{}'.format( ''.join(len(columns)*['-']) ))
        for i in mm.bonds_via_dist.nb_count:
            row = '{:<6}'.format(i)
            for count in mm.bonds_via_dist.nb_count[i]:
                row += '{:^21}'.format(mm.bonds_via_dist.nb_count[i][count])
            log.out('{}'.format(row))
     
    ####################################
    # Write molecules/cluster findings #
    ####################################
    maxID = 100 # to stop printing after maxID
    # Write molecules table          
    log.out('\n\n--------------------------------------------Cluster Analysis-------------------------------------')
    log.out('{:^10} {:^15} {:^20} {:^15} {:^15} {:^15}'.format('molID', 'Molecule Size', 'Mass', '%Mass', '%Size', 'Molecule Formula'))
    log.out('-------------------------------------------------------------------------------------------------')  
    for i in mm.molecules.data:
        data = mm.molecules.data[i]
        size = '{: >6}'.format(data.size)
        mass = '{:.2f}'.format(data.mass)
        pmass = '{:.2f}'.format(data.pmass)
        psize = '{:.2f}'.format(data.psize)
        formula = '{:^10}'.format(data.formula)
        if i <= maxID:
            log.out('{:^10} {:^15} {:^20} {:^15} {:^15} {:^15}'.format(i, size, mass, pmass, psize, formula))
    # if len(mm.molecules.data) > maxID print ...
    if len( mm.molecules.data) > maxID:
        size = '{: >6}'.format('.'); mass = '{:.2}'.format('.')
        pmass = '{:.2}'.format('.'); psize = '{:.2}'.format('.')
        formula = '{:^10}'.format('.')
        log.out('{:^10} {:^15} {:^20} {:^15} {:^15} {:^15}'.format('.', size, mass, pmass, psize, formula))
        log.out('{:^10} {:^15} {:^20} {:^15} {:^15} {:^15}'.format('.', size, mass, pmass, psize, formula))
        log.out('{:^10} {:^15} {:^20} {:^15} {:^15} {:^15}'.format('.', size, mass, pmass, psize, formula))
        
        
    ###################################################################
    # Custom extent of reaction section to write or not (hidden flag) #
    #########################################################################################################
    # Extent of reaction:                                                                                   #
    #    p = (2*(N0 - N))/(N0*fav)                                                                          #
    #                                                                                                       #
    # Critical Extent of reaction (gel point):                                                              #
    #   pg = 2/fav                                                                                          #
    #                                                                                                       #
    # Degree of polymerization:                                                                             #
    #    Xn = 2/(2 - p*fav)                                                                                 #
    #                                                                                                       #
    # Where:                                                                                                #
    #    p   = Extent of reaction                                                                           #
    #    pg  = Critical Extent of reaction (gel point)                                                      #
    #    N0  = Number of intial molecules before polymerization                                             #
    #    N   = Number of molecules at any stage or polymerization                                           #
    #    fav = Average number of functional groups present per monomer unit                                 #
    #    Xn  = Degree of polymerization                                                                     #
    #                                                                                                       #
    # Determining fav for EPON 862 example:                                                                 #
    #    1 DETDA molecule (4-functional groups)                                                             #
    #    2 DGEBF molecule (2-fuctional groups)                                                              #
    #    mix = 1*DETDA + 2*DGEBF (2 DETDA to every 1 DGEBF)                                                 #
    #    fav = (1-DETDA*4-functional-groups + 2-DGEBF*2-functional-groups)/3-molecules                      #
    #    fav = (1*4 + 2*2)/3 = 2.66                                                                         #
    #########################################################################################################
    write_extent_of_reaction = False # True to write, False to not write
    if write_extent_of_reaction:
        N = len(mm.molecules.data)
        N0 = 100    # MUST UPDATE BASED ON YOUR SYSTEM
        fav = 2.0   # MUST UPDATE BASED ON YOUR SYSTEM
        
        # Compute metrics
        p = 0; pg = 0; Xn = 0;
        if N0 > 0 and fav > 0:
            p = (2*(N0-N))/(N0*fav) # extent of reaction
            pg = 2/fav
        if p > 0 and fav > 0:
            Xn = 2/(2-p*fav) # degree of polymerization
            
        # Write metrics
        log.out('\n------------------------Curing metrics-----------------------')
        log.out('{} {:.4f}'.format('Extent of reaction aka converison      (p)   : ', p))
        log.out('{} {:.4f}'.format('Critical extent of reaction gel point  (pg)  : ', pg))
        log.out('{} {:.4f}'.format('Degree of Polymerization               (Xn)  : ', Xn))
        log.out('PLEASE NOTE: If deleting molecules using the deleta_atoms,')
        log.out('dictionary these results will be wrong!')
        
        
    #####################
    # Print By products #
    #####################
    # Write delete atoms criteria
    log.out('\n\nBy-products criteria: {}'.format(str(bp.delete_atoms)))
    
    # Warn if delete_atoms was used and files extension was not .dat or .data
    if bp.delete_atoms['criteria'] != 0 and not mm.filename.endswith('data') or mm.filename.endswith('dat'):
        log.out('WARNING Using delete_atoms option for a Non-LAMMPS file!\n')
        
    # Write by-products table
    log.out('------By Products Tally------')
    log.out('{:<14} {:>14}'.format('Type', 'Count'))
    log.out('-----------------------------')
    for i in bp.kept_molecules:
        log.out('{:<14} {:>14}'.format(i, bp.kept_molecules[i]))
        


        
        
    ###########################
    # Write ringed atoms data #
    ###########################
    # Write all inputs of find_pyramidalization
    log.out('\n\n----Inputs used for find_rings----')
    log.out('{} {}'.format('Walked along elements  : ', mm.find_rings['elements2walk']))
    log.out('{} {}'.format('Checked for ring sizes : ', mm.find_rings['rings2check']))
    log.out('{} {}'.format('Total rings found      : ', mm.rings.total))
    
    log.out(f'{mm.rings.partitioned_count} atoms along a ring seam had to be partitioned amoung the ring') 
    log.out('types (To find accurate %Mass of ring type to entire system).')
    log.out('Giving preference of partionioning in this order:')
    log.out('- 6 member ring')
    log.out('- 5 member ring')
    log.out('- 7 member ring')
    log.out('- 4 member ring')
    log.out('- 3 member ring')
    log.out('- 8 member ring')
    log.out('- minimum ring size')
    log.out('*NOTE: If count of rings exists, but no atoms exist for that ring, this means the')
    log.out('atoms for that ring were partionted to other rings that the atoms belong to.*\n')                                                                                                     
    for i in mm.rings.data:
        data = mm.rings.data[i]
        count = '{:^5}'.format(data.count)
        pcount = '{:.2f}'.format(data.pcount)
         
        # If count is greater then zero write
        if data.count > 0:
            log.out('---------------------------------------------------------------------------------------------')
            log.out('|{:^25} | {:^25} | {:^35}|'.format('Ring', 'Count', '%Ring count'))
            log.out('|{:^25} | {:^25} | {:^35}|'.format('Type', 'of Rings', 'per all rings'))
            log.out('---------------------------------------------------------------------------------------------')
            log.out('|{:^25} | {:^25} | {:^35}|'.format(i, count, pcount))
            log.out('---------------------------------------------------------------------------------------------')
            log.out('|{:^16} | {:^15} | {:^16} | {:^16} | {:^16}|'.format('Element', 'natoms', 'Mass', '%Mass', '%natoms'))
            log.out('---------------------------------------------------------------------------------------------')
            for j in mm.rings.data[i].partitioned:
                data1 = mm.rings.data[i].partitioned[j]
                size = '{:^5}'.format(data1.size)
                mass = '{:.2f}'.format(data1.mass)
                pmass = '{:.2f}'.format(data1.pmass)
                mass = '{:.2f}'.format(data1.mass)
                psize = '{:.2f}'.format(data1.psize)
                log.out('|{:^16} | {:^15} | {:^16} | {:^16} | {:^16}|'.format(j, size, mass, pmass, psize))
            log.out('---------------------------------------------------------------------------------------------\n\n')


    #################################
    # Write ringed cluster findings #
    #################################
    if len(mm.rings.clusters) > 0:
        # For new table tally based on ring type
        ring_tally = {} # { ring type : [count, ring-size, size-tally, mass-tally, pmass-tally, psize-tally] }
        for i in mm.rings.clusters:
            data = mm.rings.clusters[i]
            if data.formula in ring_tally:
                ring_tally[data.formula][0] += 1
                ring_tally[data.formula][1] = data.size
                ring_tally[data.formula][2] += data.size
                ring_tally[data.formula][3] += data.mass
                ring_tally[data.formula][4] += data.pmass
                ring_tally[data.formula][5] += data.psize
            else:
                ring_tally[data.formula] = [1, data.size, data.size, data.mass, data.pmass, data.psize]
        
        # Write new ring table
        log.out('----------------------------------Ringed Clusters------------------------------------')
        log.out('{:^10} {:>12} {:>14} {:>18} {:>16}'.format('Ring-Formula', 'Ring-Size', 'Ring-count', 'Ring-Mass-tally', 'natoms-tally'))
        log.out('-------------------------------------------------------------------------------------')  
        for formula in ring_tally:
            count, ringsize, size, mass, pmass, psize = ring_tally[formula]
            count = '{: >6}'.format(count)
            ringsize = '{: >6}'.format(ringsize)
            size = '{: >8}'.format(size)
            mass = '{:.2f}'.format(mass)
            pmass = '{:.2f}'.format(pmass)
            psize = '{:.2f}'.format(psize)
            formula = '{:<10}'.format(formula)
            log.out('{:^10} {:>12} {:>14} {:>18} {:>16}'.format(formula, ringsize, count, mass, size))

            
    #############################################
    # Write fused ring findings if user desired #
    #############################################
    if mm.find_rings['fused-rings']:
        maxID = 10 # to stop printing after maxID
        # Write ringed clusters table          
        log.out('\n\n--------------------------------------------------Fused Ring Clusters-----------------------------------------------------')
        log.out('{:<6} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>35}'.format('FusedID', 'Size', 'Mass', '%Mass', '%Size', 'Nrings', '%Rings', 'FusedRing Formula'))
        log.out('--------------------------------------------------------------------------------------------------------------------------')  
        for i in mm.rings.fused.data:
            data = mm.rings.fused.data[i]
            size = '{: >6}'.format(data.size)
            mass = '{:.2f}'.format(data.mass)
            pmass = '{:.2f}'.format(data.pmass)
            psize = '{:.2f}'.format(data.psize)
            nrings = '{:>6}'.format(data.nrings)
            prings = '{:>6}'.format(data.prings)
            formula = '{:^10}'.format(data.formula)
            if i <= maxID and data.nrings > 0:
                log.out('{:<6} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>35}'.format(i, size, mass, pmass, psize, nrings, prings, formula))
        # if len(mm.rings.clusters.data) > maxID print ...
        if len(mm.rings.fused.data) > maxID:
            size = '{: >6}'.format('.'); mass = '{:.2}'.format('.'); pmass = '{:.2}'.format('.'); psize = '{:.2}'.format('.');
            nrings = '{:>6}'.format('.');  prings = '{:>6}'.format('.'); formula = '{:^10}'.format('.');
            for i in range(3):
                log.out('{:<6} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>35}'.format('.', size, mass, pmass, psize, nrings, prings, formula))


    ###############################
    # Write hybridization results #
    ###############################
    log.out('\n\n\n-----------------------------Hybridization Information-------------------------------')
    log.out('{:^20} {:^16} {:^16} {:^16} {:^16}'.format('Atom-Type', 'natoms', 'Mass', '%Mass', '%natoms')) 
    log.out('-------------------------------------------------------------------------------------')
    for element in mm.hybridization:
        for hybridization in mm.hybridization[element]:
            atomtype = '{}-{}'.format(hybridization, element)
            size = mm.hybridization[element][hybridization].size
            mass = '{:.2f}'.format(mm.hybridization[element][hybridization].mass)
            pmass = '{:.2f}'.format(mm.hybridization[element][hybridization].pmass)
            psize = '{:.2f}'.format(mm.hybridization[element][hybridization].psize)
            
            # Only write data if natoms > 0 and do not write all hybridization
            if size > 0:# and hybridization != 'all':
                log.out('{:^20} {:^16} {:^16} {:^16} {:^16}'.format(atomtype, size, mass, pmass, psize))    

            
    #############################
    # Write atom typing results #
    #############################
    # Function to create chunks and create print string
    def divide_chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i:i + n]
    
    def pretty_print(lst):
        string = '';
        for n, i in enumerate(lst):
            # String some white space at start of string
            if n == 0: string += '{:<8}'.format('')
            string += '{:<15}'.format(i)
        return string
    
    # Print info
    log.out('\n\n\n')
    log.out('-----------------------------------------------------------------------------------------')
    log.out('| Currently supported atom types (If you know your system has an atom type that is not  |')
    log.out('| listed the force field specifc typing module needs to be added onto for that type).   |')
    log.out('| If the atom type has a trailing (T) or (F) it means that atom type has flag settings  |')
    log.out('| in the force field specifc typing module and lets you know the status of those flags. |')
    log.out('| If the atom type has a trailing (Q) it means that the written datafile has the charge |')
    log.out('| set for that atom type already (Mainly for PCFF-IFF or Metal/Mineral based FFs). For  |')
    log.out('| most atom types with the trailing (Q) there will a flag in the specific atom typing   |')
    log.out('| module to turn on or off this functionality (to give a more control of the code).     |')
    count = 0
    for i in mm.supported_types:
        types = list(mm.supported_types[i])
        count += len(types)
        if len(types) > 0:
            chunks = list(divide_chunks(types, 5))
            log.out('|-------------------------------------{:^12}--------------------------------------|'.format(i))
            for lst in chunks:
                log.out('| {:<85} |'.format(pretty_print(lst)))
    log.out('-----------------------------------------------------------------------------------------')
    log.out('| {:^85} |'.format(f'Total count of supported atom types: {count}'))
    log.out('-----------------------------------------------------------------------------------------')
    
    ###############################
    # print results to the screen #
    ###############################
    type_percent = '{:.2f}'.format(100*mm.tally["found"]/len(mm.atoms))
    log.out('\n\n\nFinal outcome of found atom types:')
    log.out(f'Total atom count                   : {len(mm.atoms)}')
    log.out(f'Total Parameterized atom count     : {mm.tally["found"]}')
    log.out(f'Assumed Parameterized atom count   : {mm.tally["assumed"]}')
    log.out(f'Failed Parameterized atom count    : {mm.tally["failed"]}')
    log.out(f'Total Parameterized atom perentage : {type_percent}')
    
    
    ###################################################
    # Find count of each atom type assigned and write #
    ###################################################
    atomtypes_lst = sorted({mm.atoms[i].nta_type for i in mm.atoms})
    atomtypes_count = {i:0 for i in atomtypes_lst}
    for i in mm.atoms:
        atomtypes_count[mm.atoms[i].nta_type] += 1

    # Write atom types table
    log.out('\n------Atom types Tally------')
    log.out('{:<16}  {:>10}'.format('Type', 'Count'))
    log.out('----------------------------')
    for i in atomtypes_count:
        log.out('{:<16}  {:>10}'.format(i, atomtypes_count[i]))
        
    return