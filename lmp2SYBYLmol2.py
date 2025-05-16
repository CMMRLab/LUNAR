# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
May 16, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

Will read in a LAMMPS datafile to convert to a SYBYL
.mol2 file. This is useful for "drawing" in atoms 
and bonds into a LAMMPS datafile. Avogadro2 seems
to be the best software to read in the .mol2 files
this code writes. ChemDraw can also read the .mol2
file, if there are fewer atoms (~<1000). The written
.mol2 file can also be opened in VMD as well.

    **************************************************************
    * Requirements:                                              *
    *   python 3.7+                                              *
    *                                                            *
    * Run methods:                                               *
    *   - IDE (manipulate variables and run from IDE)            *
    *   - GUI (manipulate variables and run from. Default        *
    *          settings set from this script)                    *
    *   - command line (python3 lmp2SYBYLmol2.py -man to get     *
    *                   command line override options and        *
    *                   examples)                                *
    *                                                            *
    *                                                            *
    **************************************************************
"""


##############
### Inputs ###
##############
##################################################################################################################
# Python boolean variable to use GUI or not (True or False). All GUI settings are intialized from this script so #
# adjusting this script will set the default settings that the GUI will load with. Please NOTE that when using   #
# the GUI to check the console or terminal print outs every time a system is run through the code because that is#
#  where the import information will be displayed. Examples:                                                     #
#   use_GUI = True  # Will launch GUI to manipulate variables and run from GUI                                   #
#   use_GUI = False # Will NOT launch GUI where to manipulate variables and run will have to be from the IDE     #
#                                                                                                                #
# Python int value to range from 0 to 200 to control the GUI zoom. Where it option was added to account for the  #
# different ways individuals have their OS display settings setup. A GUI_zoom = 100, means that default settings #
# are used, whereas a GUI_zoom = 120 means increase GUI size by 20% and GUI_zoom = 80 means decrease GUI by 20%. #
# Examples:                                                                                                      #
#   GUI_zoom = 100 # use default GUI size                                                                        #
#   GUI_zoom = 75  # decrease GUI size by 25%                                                                    #
#                                                                                                                #
# Update use_GUI and GUI_zoom as desired.                                                                        #
##################################################################################################################
use_GUI = True
GUI_zoom = 100


##################################################################################################################
# LAMMPS topofile with bonds in it (THIS IS A MUST!). The LAMMPS datafile can have new LAMMPS "type labels" in   #
# it. The only supported atom styles are full, charge, or molecular. More atom styles can be coded if needed.    #
##################################################################################################################
topofile ='EXAMPLES/lmp2SYBYLmol2/detda_typed_IFF.data'


##################################################################################################################
# Python string variable type to set parent directory to store all new files this code will write. If the        #
# variable is left as an empty string, files will be written to the present working directory. Example to set    #
# parent_directory as present working directory:                                                                 #
#     parent_directory = '' or parent_directory = '.'                                                            #
#                                                                                                                #
# Example to set parent_directory as with path from topofile:                                                    #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile'                               #
#     files will be written to 'Furan_Resin/files/' since the 'topofile' string envokes this option.             #
#                                                                                                                #
# Example to set parent_directory as with path from topofile and build dirs from that location:                  #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile/NEWDIR'                        #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/NEWDIR' will build relative        #
#     directories                                                                                                #
#                                                                                                                #
#     topofile = 'Furan_Resin/files/furan.data' and  parent_directory = 'topofile/../NEWDIR'                     #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/../NEWDIR' will build relative     #
#     directories                                                                                                #
#                                                                                                                #
# Example to set parent_directory to a location inside of present working dirctory:                              #
#     parent_directory = 'inside'       # will save files in $pwd/inside                                         #
#     parent_directory = 'inside/test'  # will save files in $pwd/inside/test                                    #
#                                                                                                                #
# Example to set parent_directory to a location outside of present working dirctory:                             #
#     parent_directory = '../'      # will save files in $pwd/../                                                #
#     parent_directory = '../test'  # will save files in $pwd/../test                                            #
#                                                                                                                #
# Update parent_directory as desired.                                                                            #
##################################################################################################################
parent_directory = 'topofile'


##################################################################################################################
# Option to remove periodic bonds from the written .mol2 file. Typically you would only use this option to       #
# visualize the converted files in ChemDraw, VMD, Avogadro, Avogadro2, ... since periodic bonds will be          #
# visualized as spanning the imaginary simulation cell. If you are using lmp2SYBYLmol2 for adding atoms for more #
# MD simulations DO NOT USE this option (True or False).                                                         #
##################################################################################################################
remove_PBC_bonds = False


##################################################################################################################
# A Boolean to add a "pseudo simulation cell" to the written .mol2 file. The "pseudo simulation cell" is added   #
# to the file by setting the box corners as atoms and the box edges as bonds. The purpose of this option is to   #
# allow for the visualization of your model in VMD with a simulation cell defined, since VMD does not natively   #
# support simulation cells nor does the .mol2 file format. When the atoms are written to the .mol2 file the      #
# following attributes are set:                                                                                  #
#   - column1 -> atomID:      incremented atomID                                                                 #
#   - column2 -> element:     set as "Bx" (shorthand for Box)                                                    #
#   - column3 -> x:           box corner X-postion                                                               #
#   - column4 -> y:           box corner Y-postion                                                               #
#   - column5 -> z:           box corner Z-postion                                                               #
#   - column6 -> atom type:   set as "Bx" (shorthand for Box)                                                    #
#   - column7 -> subst_id:    set as max(molIDs)+1 (box atoms will always have a molID greater than any atom)    #
#   - column8 -> subst_name:  set as "BOX"                                                                       #
#   - column9 -> charge:      set as "0.000" as charge is meaningless here                                       #
#                                                                                                                #
# Knowing which attributes are set and how they are set allows you to generate different visualizations in VMD.  #
# For example the attribute names map onto VMD coloring and selection categories as such:                        #
#   subst_id   -> resid                                                                                          #
#   subst_name -> resname                                                                                        #
#                                                                                                                #
# As an example, set addbox=True for the default "EXAMPLES/lmp2SYBYLmol2/detda_typed_IFF.data" topofile, run,    #
# and look at the bottom of the @<TRIPOS>ATOM and @<TRIPOS>BOND section, to see the "box" atoms and bonds. Then  #
# do the following in VMD:                                                                                       #
#   1. File -> New Molecule -> Browse -> Select "detda_typed_IFF.mol2" file -> Load                              #
#   2. Graphics -> Representations  &&  Graphics -> Colors                                                       #
#     * Which will open both the Representations and Colors dialog boxes                                         #
#      * In the Representations dialog box do the following:                                                     #
#        * Set Selected atoms as "name C H N", which generates a representation of just the molecule             #
#        * Change the Coloring method to "Type"                                                                  #
#        * Change the Drawing method to "CPK"                                                                    #
#                                                                                                                #
#        * Click "Create rep" to generate another representation Do the following for this representation        #
#          * Set Selected atoms as "resname BOX", which generates a representation of just the atoms of the box  #
#          * Change the Coloring method to "ResName"                                                             #
#          * Change the Drawing method to "Lines"                                                                #
#                                                                                                                #
#       * In the Colors dialog box do the following:                                                             #
#         * Categories -> Type  && Names -> C  && Colors -> black                                                #
#         * which changes the carbon atoms to color black. Can color the remaining atoms in the same way         #
#         * Categories -> Resname  && Names -> BOX  && Colors -> blue                                            #
#           * which changes the box atoms to the color blue                                                      #
#           * NOTE: VMD has predefined Resname's and any new resnames that are defined when reading a file are   #
#                   set at the bottom of the list, so you can just assume this and scroll to the bottom of the   #
#                   Names menu, to find "BOX"                                                                    #
# 	                                                                                                             #
# There are a few other ways to generate your representation's, however those are self studies. A few useful     #
# entries for the "Selected Atoms" box are (in example format, change numbers and letters accoridingly):         #
#   all and x > -12 and y > -12 and z > -12                                                                      #
#   resid > 10                                                                                                   #
#   name C O                                                                                                     #
#   mass 5 to 11.5                                                                                               #
#   index < 10                                                                                                   #
#   within 5 of name H                                                                                           #
##################################################################################################################
addbox = False


##################################################################################################################
# The mass_map dictionary is now a "global" dictionary stored in src/masses.py. The purpose of this was to       #
# simplify adding new elements, where the new elements can now be applied to every code that uses the mass_map.  #
# If you get an "ERROR Not all masses in ... are in the mass_map dictionary.", you will now have to open         #
# src/masses.py and update the mass_map dictionary found in that file.                                           #
##################################################################################################################
import src.masses as masses
mass_map = masses.mass_map 




#################################
### Main conversion function ####
#################################
def main(topofile, parent_directory, remove_PBC_bonds, mass_map, addbox, log=None):
    ##############################
    # Import Necessary Libraries #
    ##############################
    import src.read_lmp as read_lmp
    import src.io_functions as io_functions
    from collections import OrderedDict
    import sys
    import os
    
    # Configure log (default is level='production', switch to 'debug' if debuging)
    if log is None:
        log = io_functions.LUNAR_logger()
    log.configure(level='production')
    #log.configure(level='debug')
    
    
    ########################################################
    # set version and print starting information to screen #
    ########################################################
    version = 'v1.2 / 16 May 2025'
    log.out(f'\n\nRunning lmp2SYBYLmol2 {version}')
    log.out(f'Using Python version {sys.version}')
    
    
    ##########################################################
    # Read LAMMPS datafile into memory as class "m" applying #
    # forward mapping of type labels (if applicable)         #
    ##########################################################
    if os.path.isfile(topofile):
        m = read_lmp.Molecule_File(topofile, method='forward', sections=['Atoms', 'Bonds'])
        log.out(f'Read in {m.filename} LAMMPS datafile')
    else: log.error(f'ERROR lammps datafile: {topofile} does not exist')
    basename = os.path.basename(topofile) # Strip path from topofile
    if basename.endswith('data.gz'):
        basename = basename.replace('.gz', '') # remove .gz suffix
    basename = basename[:basename.rfind('.')] # Find file basename
    log.debug(f'basename = {basename}')
    
    
    ##################################
    # Function to get element symbol #
    ##################################
    def get_element_symbol(m, atomID, mass_map):
        mass = m.masses[m.atoms[atomID].type].coeffs[0]
        # Try gettint element symbol
        element = 'DEFAULT' # Initialize and update if found
        try: element = [i for i in mass_map if mass in mass_map[i]][0]
        except: log.error(f'ERROR Not all masses in {topofile} are in the mass_map dictionary. Failed for mass: {mass}')
        return element
    
    
    #######################################
    # Setting up directories and where to #
    # write final files and results to    #
    #######################################
    # Find present working directory and find/create paths to store code results
    pwd = os.getcwd()
    path = os.path.join(pwd, parent_directory)
    
    # If parent_directory == 'topofile' use topofile path as path
    if 'topofile' in parent_directory:
        log.out('Using path from topofile to set parent_directory ...')
        path = io_functions.get_dir_from_topofile(topofile, parent_directory)
    
    # Check if path exits. IF not create
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)
        
    # Change the current working directory to path so all files get written to that directory
    os.chdir(path)
    
    
    ################################################################
    # Check if bond is periodic using the minimum image convention #
    ################################################################
    # Find box dimensions to remove periodic boundary conditions
    x = m.xbox_line.split(); y = m.ybox_line.split(); z = m.zbox_line.split();
    xlo = x[0]; xhi = x[1]; ylo = y[0]; yhi = y[1]; zlo = z[0]; zhi = z[1];
    lx = float(x[1])-float(x[0]); ly = float(y[1])-float(y[0]); lz = float(z[1])-float(z[0]);
    
    # set max_x, max_y, max_z w/ minimum image convention
    max_x = lx/2; max_y = ly/2; max_z = lz/2;
    
    # Function to check bond periodicity status
    def check_bond_periodicity(m, id1, id2):
        pbc_flag = False # Intialize and update if bond is periodic
        x1 = m.atoms[id1].x; y1 = m.atoms[id1].y; z1 = m.atoms[id1].z
        x2 = m.atoms[id2].x; y2 = m.atoms[id2].y; z2 = m.atoms[id2].z
        
        # if bond fails minimum image convention it is periodic
        if abs(x2 - x1) > max_x: pbc_flag = True
        if abs(y2 - y1) > max_y: pbc_flag = True
        if abs(z2 - z1) > max_z: pbc_flag = True
        return pbc_flag
    
    # Find bonds that we want to write
    bondIDs2write = []
    for i in m.bonds:
        id1, id2 = m.bonds[i].atomids
        if remove_PBC_bonds: pbc_flag = check_bond_periodicity(m, id1, id2)
        else: pbc_flag = False
        if not pbc_flag: bondIDs2write.append(i)
        
        
    # atom positons and bonds to define simulation cell for VMD
    boxatoms = [(xhi, ylo, zhi), (xlo, ylo, zhi), (xlo, yhi, zhi), (xhi, yhi, zhi),
                (xhi, ylo, zlo), (xlo, ylo, zlo), (xlo, yhi, zlo), (xhi, yhi, zlo)]
    boxbonds = [(1, 2), (2, 3), (3, 4), (4, 1), (5, 6), (6, 7),
                (7, 8), (8, 5), (4, 8), (5, 1), (3, 7), (6, 2)]


    ##########################################################################################
    # Writing new file mol2 with bonds information                                           #
    # https://chemicbook.com/2021/02/20/mol2-file-format-explained-for-beginners-part-2.html #
    ##########################################################################################
    with open(basename+'.mol2','w') as f: 
        # Write molecule section
        f.write('@<TRIPOS>MOLECULE\n')
        f.write(f'{m.header} > lmp2SYBYLmol2 {version} w/remove_PBC_bonds={str(remove_PBC_bonds)}\n')
        if not addbox:
            f.write(f'  {len(m.atoms)} {len(bondIDs2write)}    0    0    0\n')
        else:
            f.write(f'  {len(m.atoms)+len(boxatoms)} {len(bondIDs2write)+len(boxbonds)}    0    0    0\n')
        f.write('SMALL\n')
        f.write('NO_CHARGES\n')
        f.write('****\n')
        f.write('Energy = 0\n')
        
        # Write Atoms info
        f.write('\n@<TRIPOS>ATOM\n')
        m.atoms = dict(OrderedDict(sorted(m.atoms.items()))) # sort to keep IDs as close as possible to orginal
        id_map = {} # { orginal atomID : New atomID } to make IDs contiguous if not already
        molids = set([1])
        for n, i in enumerate(m.atoms, 1):
            atom = m.atoms[i]
            
            # Find atoms info
            element = get_element_symbol(m, i, mass_map)
            x = '{:>17.4f}'.format(atom.x) # float point
            y = '{:>10.4f}'.format(atom.y) # float point
            z = '{:>10.4f}'.format(atom.z) # float point
            try:
                subst_id = atom.molid # The ID number of the substructure containing the atom (int)  [VMD RESID Coloring]
                molid = atom.molid
            except:
                subst_id = 1 # The ID number of the substructure containing the atom (int)  [VMD RESID Coloring]
                molid = 1
            molids.add(molid)
            subst_name = '****' # The name of the substructure containing the atom (string)
            subst_name = '{:>4}'.format('m'+str(molid)) # Will give access through [VMD ResName Coloring]
            charge = '{:>10.4f}'.format(atom.charge)
            
            # Add id to map
            id_map[i] = n
            
            # Write atoms info
            f.write('{:>7} {:<2} {} {} {} {:<2} {:>7} {:>7} {}\n'.format(n, element, x, y, z, element, subst_id, subst_name, charge))
            
        # Add in box atoms if flag
        if addbox:
            boxID = max(molids) + 1; box_map = {} # { index of location : new atomID }
            for ID, box in enumerate(boxatoms, 1):
                n += 1 
                x, y, z = box
                element = 'Bx'
                box_map[ID] = n
                x = '{:>17.4f}'.format(float(x)); y = '{:>10.4f}'.format(float(y))
                z = '{:>10.4f}'.format(float(z)); charge = '{:>10.4f}'.format(0);
                subst_name = 'BOX' 
                subst_id = boxID
                f.write('{:>7} {:<2} {} {} {} {:<2} {:>7} {:>7} {}\n'.format(n, element, x, y, z, element, subst_id, subst_name, charge))
            
        # Write Bonds info
        f.write('@<TRIPOS>BOND\n')
        bondIDs2write = sorted(bondIDs2write) # sort to keep IDs as close as possible to orginal
        for n, i in enumerate(bondIDs2write, 1):
            id1, id2 = m.bonds[i].atomids
                        
            # Set bond_type as 1 for now...
            # Possible options:
            #    1 = single
            #    2 = double
            #    3 = triple
            #    am = amide
            #    ar = aromatic
            #    du = dummy
            #    un = unknown (cannot be determined from the parameter tables)
            #    nc = not connected
            bond_type = '1'
            
            # Write bonds info
            new_id1 = id_map[id1]; new_id2 = id_map[id2]
            f.write('{:>6} {:>6} {:>6} {:>6}\n'.format(n, new_id1, new_id2, bond_type))
            
        if addbox:
            for id1, id2 in boxbonds:
                new_id1 = box_map[id1]
                new_id2 = box_map[id2]
                bond_type = 'du'
                n += 1;
                f.write('{:>6} {:>6} {:>6} {:>6}\n'.format(n, new_id1, new_id2, bond_type)) 
            
    #######################
    # Clean up and ending #
    #######################
    # Print file locations and completion of code
    log.out(f'\n\nAll outputs can be found in {path} directory')
    log.out('\n\nNormal program termination\n\n')
    
    # Show number of warnings and errors
    log.out_warnings_and_errors()
    
    # write log
    log.write_logged(basename+'.log.lunar')
    
    # Change back to the intial directory to keep directory free for deletion
    os.chdir(pwd)
    return

###################################
### Import needed files and run ###
###################################
if __name__ == "__main__":  
    import src.command_line as cl
    import sys

    ####################################################################################################
    # Set IDE/command line mode option:                                                                #
    # - if run in IDE commandline_inputs lst will have zero entries. Which means use hard coded inputs #
    #   from inputs section.                                                                           #
    #                                                                                                  #
    # - if run in command line, but commandline_inputs still has zero entries, this means use hard     #
    #   coded inputs from inputs section.                                                              #
    #                                                                                                  #
    # - else run in command line and commandline_inputs has entries determine what entries will be     #
    #   overrided from hard coded inputs from inputs section. This overide will occur in in the inputs #
    #   section called Command Line Override. *NOTE: if user does not specify certain command line     #
    #   options, the default settings from inputs section will be enforced. A message will tell user   #
    #   if inputs from inputs section have been enforced. The user can find all command line options   #
    #   by running:                                                                                    #
    #       python3 lmp2SYBYLmol2.py -opt                                                              #
    #              or                                                                                  #
    #       python3 lmp2SYBYLmol2.py -man                                                              #
    #   in the terminal and the code will print out all command line option and terminate before any   #
    #   further analysis is done                                                                       #
    ####################################################################################################
    # Find command line inputs that are not the python script itself
    commandline_inputs = [i for n, i in enumerate(sys.argv) if n > 0]
    if commandline_inputs and '-gui' not in commandline_inputs: 
        topofile, parent_directory, remove_PBC_bonds = cl.lmp2SYBYLmol2_interface(topofile, parent_directory, remove_PBC_bonds, commandline_inputs)
        use_GUI = False

    if use_GUI or '-gui' in commandline_inputs:
        print('\n\n\nlmp2SYBYLmol2 is currently running in GUI mode, where all GUI inputs are intialized from lmp2SYBYLmol2.\n\n\n')
        from src.lmp2SYBYLmol2_GUI import GUI 
        GUI(topofile, parent_directory, remove_PBC_bonds, mass_map, addbox, GUI_zoom)
    else: main(topofile, parent_directory, remove_PBC_bonds, mass_map, addbox)
