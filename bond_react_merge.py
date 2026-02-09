# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.16
November 2nd, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

    **********************************************************
    * Requirements:                                          *
    *   python 3.7+                                          *
    *                                                        *
    * Dependencies:                                          *
    *   datafiles were generated with all2lmp or were        *
    *   commented by bond_react_merge_prep.py                *
    *                                                        *
    * Run methods:                                           *
    *   - IDE (manipulate variables and run from IDE)        *
    *   - GUI (manipulate variables and run from. Default    *
    *          settings set from this script)                *
    *   - command line (python3 bond_react_merge.py -man to  *
    *                   get command line override options and*
    *                   examples)                            *
    *                                                        *
    **********************************************************
"""


##############
### Inputs ###
##############
#############################################################################################################################
# Python boolean variable to use GUI or not (True or False). All GUI settings are intialized from this script so adjusting  #
# this script will set the default settings that the GUI will load with. Please NOTE that when using the GUI to check the   #
# console or terminal print outs every time a system is run through the code because that is where the import information   #
# will be displayed. Examples:                                                                                              #
#   use_GUI = True  # Will launch GUI to manipulate variables and run from GUI                                              #
#   use_GUI = False # Will NOT launch GUI where to manipulate variables and run will have to be from the IDE or command line#
#                                                                                                                           #
# Python variable to set the maxfiles the GUI will load with (maxfiles must be an integer value). Please NOTE that the GUI  #
# supports file overloading when running, so maxfiles only sets the intial number of files and more will be dynamically     #
# added if more files are provide while the GUI is loaded. Examples:                                                        #
#   maxfiles = 4 # GUI will load with 4 files slots and the GUI will automatically add more when "overloaded"               #
#   maxfiles = 8 # GUI will load with 8 files slots and the GUI will automatically add more when "overloaded"               #
#                                                                                                                           #
# Python int value to range from 0 to 200 to control the GUI zoom. Where it option was added to account for the different   #
# ways individuals have their OS display settings setup. A GUI_zoom = 100, means that default settings are used, whereas a  #
# GUI_zoom = 120 means increase GUI size by 20% and GUI_zoom = 80 means decrease GUI by 20%. Examples:                      #
#   GUI_zoom = 100 # use default GUI size                                                                                   #
#   GUI_zoom = 75  # decrease GUI size by 25%                                                                               #
#                                                                                                                           #
# Update use_GUI, GUI_zoom, scroll_bar, and maxfiles as desired.                                                            #
#############################################################################################################################
use_GUI = True
GUI_zoom = 100
maxfiles = 6


#############################################################################################################################
# The main input to this code is a python dictionary where the keys represent a "file-tag". The key (file-tag) will be what #
# is used to decide what type of final file to write and what type of processing needs to be done. The value of the         #
# dictionary will be a string for calling the LAMMPS datafile you want to read. Rules for formating the files dictionary:   #
# - the key can be any of the following (depending on desired output):                                                      #
#    * data1, data2, ... dataN for an arbitrary amount of starting molecules                                                #
#    * pre1, pre2, ...   preN for an arbitrary amount of pre rxns                                                           #
#    * post1, post2, ... postN for an arbitrary amount of post rxns                                                         #
#    * infile            for reading in input script that specifies file names                                              #
#    * NOTE: if generate_map_file the N-ID assigned between preN/postN MUST BE CONSISTENT!!!!                               #
# - the value MUST a LAMMPS datafile generated by all2lmp or commented by bond_react_merge_prep.py. Alternatively if the    #
#   datafile has type labels generated by all2lmp and has been run through LAMMPS already, the all2lmp style comments will  #
#   by removed, but the src/read_lmp.py code will re-insert the comments based on the type labels that LAMMPS preserves.    #
#   This will make for an easier pre/post/medial processing of merging parameters and generating templates if needed.       #
#   Please ensure the following:                                                                                            #
#    * datafile mapping to the "file-tag" is proper (otherwise files will be written incorrectly)                           #
#    * that the files keys are all unqiue (IE you dont use data1 twice) because python dictionaries can only have a set of  #
#      unqiue keys. Hence the N-ith notation.                                                                               #
#                                                                                                                           #
# Format of the files dictionary:                                                                                           #
#    files = { "file-tag" :   'name-of-tagged-datafile' }                                                                   #
#                                                                                                                           #
# Example of reading in two dataN tagged files and a single pre/post reaction file pairing (generate_map_file will find     #
# equivalences between the 'pre1' tag and the 'post1' tagged files. map_near_edge_rxn_charges will search for charges       #
# of the 'pre1' and 'post1' tagged files by comparing to all dataN tagged files):                                           #
#    files = {'data1' : detda_typed_IFF.data,     # Contains a single molecule                                              #
#             'data2' : detda_typed_IFF.data,     # Contains a single molecule                                              #
#             'pre1'  : pre_rxn1_typed_IFF.data,  # pre-reacted bond/react file (tagged: 1 )                                #
#             'post1' : post_rxn1_typed_IFF.data} # post-reacted bond/react file (tagged: 1)                                #
#                                                                                                                           #
# Example of reading in a set of files via an input script (read from $pwd/test directory):                                 #
#    files = {'infile' : 'test/example_merge_files.txt'}, where the format for an input script is as follows:               #
#             - "#" is a comment character and anything trailing after "#" will be ignored                                  #
#             - currently supported options to add in an "infile":                                                          #
#                 path="dir/layout/paths" # path to append to each file name specified (" characters must bound path)       #
#                 parent_directory="dir"  # parent directory to write files to (overrides parent directory python variable) # 
#                 dataN example1.data     # the dataN "file-tag" followed by white space and the tagged datafile            #
#                 preN  pre_rxn1.data     # the preN  "file-tag" followed by white space and the tagged datafile            #
#                 postN pre_rxn1.data     # the postN "file-tag" followed by white space and the tagged datafile            #
#                                                                                                                           #
# Update files dictionary as needed.                                                                                        #
#############################################################################################################################
files = {'data1'  : 'EXAMPLES/EPON_862/all2lmp_Outputs/detda_typed_IFF.data',           # monomer
         'data2'  : 'EXAMPLES/EPON_862/all2lmp_Outputs/dgebf_typed_IFF.data',           # monomer
         'pre1'   : 'EXAMPLES/EPON_862/all2lmp_Outputs/pre_reaction_1_typed_IFF.data',  # for rxn1
         'post1'  : 'EXAMPLES/EPON_862/all2lmp_Outputs/post_reaction_1_typed_IFF.data', # for rxn1
         'pre2'   : 'EXAMPLES/EPON_862/all2lmp_Outputs/pre_reaction_2_typed_IFF.data',  # for rxn2
         'post2'  : 'EXAMPLES/EPON_862/all2lmp_Outputs/post_reaction_2_typed_IFF.data'} # for rxn2


# Comment/uncomment and change testN dir to run test cases. *NOTE: parent_directory set in each merge_files.txt
# to help keep files in correct locations and minimize the amount of effort required to run test cases.*
#files = {'infile'  : 'EXAMPLES/bond_react_merge_tests/test11/merge_files_tag=infile.txt'}


#############################################################################################################################
# Is a Python string variable in which to set the new output filename(s). The following options exist for using the newfile #
# string for setting the output file basenames:                                                                             #
#                                                                                                                           #
#   if newfile starts with ':' or ends with ':'                                                                             #
#     The output filename(s) will be the same as the input filename(s), but will have a suffix or prefix added to the file  #
#     basename. The following are examples:                                                                                 #
#       Suffix (newfile = ':_merged'  and topofile = 'detda.data')                                                          #
#         basename = 'detda_merged', where the ':' character acts as a placeholder for the topofile basename.               #
#       Prefix (newfile = 'merged-:'  and topofile = 'detda.data')                                                          #
#         basename = 'merged-detda', where the ':' character acts as a placeholder for the topofile basename.               #
#     Recommended usage: common and safe as this method is safe and output filename(s) carry similar names to input         #
#     filename(s).                                                                                                          #
#                                                                                                                           #
#   if newfile == 'ANYTEXT'                                                                                                 #
#     The output filename(s) will be set as the file tag of each file and 'ANYTEXT' will be a suffix appended to the file   #
#     tag. For example:                                                                                                     #
#       newfile = '-merged', topofile = 'detda.data', and tag = 'data1')                                                    #
#         basename = 'data1-merged'                                                                                         #
#     Recommended usage: occasional and safe as output filename(s) no longer carry similar names as output filename(s), but #
#     is safe as output filename(s) will not overwrite input filename(s)                                                    #
#                                                                                                                           #
#   if newfile == ''                                                                                                        #
#     The output filename(s) will be the same as the input filename(s). This can be a dangerous option as you may           #
#     inadvertently overwrite a file and then must assume the file to contain certain information, but it contains other    #
#     information.                                                                                                          #
#     Recommended usage: rare and dangerous as this could lead to using incorrect files if not paying attention very        #
#     carefully.                                                                                                            #
#                                                                                                                           #
# Update files newfile as needed.                                                                                           #
#############################################################################################################################
newfile = ':_merged'


#############################################################################################################################
# Python string variable type to set parent directory to store all new files this code will write. If the variable is left  #
# as an empty string, files will be written to the present working directory. Examples:                                     #
# Example to set parent_directory as present working directory:                                                             #
#     parent_directory = '' or parent_directory = '.'                                                                       #
#                                                                                                                           #
# Example to set parent_directory as with path from 1st topofile in the files dictionary:                                   #
#     files = {'data1':'Furan_Resin/files/furan.data'} and  parent_directory = 'topofile'                                   #
#     files will be written to 'Furan_Resin/files/' since the 'topofile' string envokes this option.                        #
#                                                                                                                           #
# Example to set parent_directory as with path from 1st topofile in the files dictionary and build dirs from that location: #
#     files = {'data1':'Furan_Resin/files/furan.data'} and  parent_directory = 'topofile/NEWDIR'                            #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/NEWDIR' will build relative directories       #
#                                                                                                                           #
#     files = {'data1':'Furan_Resin/files/furan.data'} and  parent_directory = 'topofile/../NEWDIR'                         #
#     files will be written to 'Furan_Resin/files/NEWDIR' since the 'topofile/../NEWDIR' will build relative directories    #
#                                                                                                                           #
# Example to set parent_directory to a location inside of present working dirctory:                                         #
#     parent_directory = 'inside'       # will save files in $pwd/inside                                                    #
#     parent_directory = 'inside/test'  # will save files in $pwd/inside/test                                               #
#                                                                                                                           #
# Example to set parent_directory to a location outside of present working dirctory:                                        #
#     parent_directory = '../'      # will save files in $pwd/../                                                           #
#     parent_directory = '../test'  # will save files in $pwd/../test                                                       #
#                                                                                                                           #
# Update parent_directory as desired.                                                                                       #
#############################################################################################################################
#parent_directory = 'topofile/../bond_react_merge_Outputs'
parent_directory = 'topofile'


#############################################################################################################################
# Python string variable to set atom style format of the written data file. Currently supported atom styles:  charge,       #
# molecular, full, angle, bond, atomic, dipole, dpd, or line. Example:                                                      #
#     atom_style = 'full' # will set atom style in written data file as full                                                #
#                                                                                                                           #
# Update atom_style as desired.                                                                                             #
#############################################################################################################################
atom_style = 'full'


#############################################################################################################################
# Python boolean variable to generate the map file based on a cost-fitting method (True or False). Atoms in the file MUST   #
# have a comments as :  atom-type/element for this option to work (IE cp/C, h/H, op/O, ...). Examples:                      #
#    generate_map_file = True  # Will attempt to generate map file of preN/postN tagged pairs                               #
#    generate_map_file = False # Will not attempt to generate map file of preN/postN tagged pairs                           #
#                                                                                                                           #
# Additional information for using generate_map_file option can be specified in the preN tagged file header for the Nth-rxn #
# NOTE that the code will only look for this information in the preN tagged file! Some are optional and others are required #
# depending on the reaction trying to be modeled. The following are supported additional options to put in the preN tagged  #
# file header:                                                                                                              #
#    BondingIDs = [1, 2, 3, 4]                                                                                              #
#      A string that looks like a python variable and a python list. The variable name BondingIDs specifies what this info  #
#      is and must be present for usage of this option! This option is to specify "linking" atomIDs in the pre-rxn template #
#      and the post-rxn template to help the code set equivalances. NOTE that the usage of this option is NOT required and  #
#      the code will still attempt to find the correct mapping without this being specified in the preN tagged file.        #
#      However if BondingIDs are specified it will help ensure the code can find the mapping between the pre and post-rxn   #
#      by using them as starting points for a BFS (Breadth First-Search) traversal, mapping the pre-atomIDs to the          #
#      post-atomIDs at each sequential depth stage emanating from the BondingIDs. Not all rxn templates can use this option #
#      since it assumes minimal topological changes from the pre to the post-rxn, where the "linking" pair can be seen in   #
#      both the pre and post-rxn template. The BondingIDs option is equivalent to that of AutoMapper's --ba option (if you  #
#      are familar), where 4 atomIDs MUST BE specified and the order is important. The first two IDs are the IDs for the    #
#      "linking" atomIDs in the pre-rxn and the second two IDs are for the "linking" atomIDs in the post-rxn template. The  #
#      1st atomID and the 3rd atomID MUST BE equivalences between the pre and post-rxn template and the 2nd atomID and 4th  #
#      atomID MUST BE equivalences between the pre and post-rxn template. To use this option simply put the string of       #
#      BondingIDs = [1, 2, 3, 4], in the header of the corresponding pre-rxn datafile before running bond_react_merge.py.   #
#      NOTE the location in the header is not important as long as the string is in the format shown above. Examples:       #
#          HEADER > atom_typing v1.0 / 16 Jan 2023 > all2lmp: v1.4 / 9 Jan 2023 > BondingID = [7, 14, 7, 9]                 #
#          HEADER, BondingIDs = [1, 2, 3, 4] postioning does not matter as long as string format is correct                 #
#          BondingIDs = [1, 2, 3, 4] postioning does not matter as long as string format is correct                         #
#                                                                                                                           #
#    CreateIDs = [31, 32, 33]                                                                                               #
#      A string that looks like a python variable and a python list. The variable name CreateIDs specifies what this info is#
#      and must be present for usage of this option! This option is to specify the atomIDs that are considered to be atoms  #
#      that bond/react can create. These atomIDs are the atomIDs found in the post-rxn template and DO NOT have a map to    #
#      any pre-rxn atomIDs. If the rxn-template has any atomIDs that are being created, THIS MUST APPEAR in the preN tagged #
#      file for that rxn if you want to use the generate_map_file option. If the rxn creates new atomIDs and they are not   #
#      specified by this option and you try to use the generate_map_file option the code will exit with an ERROR saying the #
#      pre-rxn template and post-rxn template DO NOT have the same number of atoms. To use this option simply put the string#
#      of CreateIDs = [1, 2, 3], in the header of the corresponding pre-rxn datafile before running bond_react_merge.py.    #
#      NOTE the location in the header is not important as long as the string is in the format shown above. Examples:       #
#          HEADER > atom_typing v1.0 / 16 Jan 2023 > all2lmp: v1.4 / 9 Jan 2023 > CreateIDs = [31, 32, 33]                  #
#          HEADER, CreateIDs = [31, 32, 33]  postioning does not matter as long as string format is correct                 #
#          CreateIDs = [31, 32, 33]  postioning does not matter as long as string format is correct                         #
#                                                                                                                           #
#    Reduce = [1, 2, 3, 4, 3]                                                                                               #
#      A string that looks like a python variable and a python list. The variable name Reduce specifies what this info is   #
#      and must be present for usage of this option! If the length of the reduce list is 5 elements long with the following #
#      meaning of each index:  [preID1, preID2, postID1, postID2, depth], where preID1 and postID1 must be equivs and       #
#      preID2 and postID2 must be equivs and depth is the Breadth-First Traversal (DFS) depth of atoms to keep emanating    #
#      outward from the ID1 and ID2 locations. Usage of this option will keep all reduce down the number of atoms in the    #
#      preN/postN tagged files, where all atoms in the depth range from ID1 and ID2 are kept. The location in the header    #
#      does not matter for usage of this option. NOTE that ID1 and ID2 need not be "linking" atomIDs but should be near the #
#      "center" of the post-rxn template. Using Reduce with list length of imposes a depth restriction from ID1 and ID2, to #
#      be identical. See coupling of Reduce and BondingIDs to get around the depth restriction.                             #
#                                                                                                                           #
#    Reduce = [3, 5]; BondingIDs = [12, 33, 12, 15];  or  BondingIDs = [12, 33, 12, 15] Reduce = [3, 5] etc ...             #
#      A string that looks like a python variable and a python list. The coupling of Reduce and BondingIDs specifies what   #
#      this info is and must be present for usage of this option! If the length or reduce list is 2 elements and length of  #
#      BondingIDs is 4 elements the BondingIDs are used as the as the points to emanate out from and the 2 elements of the  #
#      reduce list set the max depth to search from for atoms to keep. When coupling both Reduce and BondingIDs options the #
#      indexes has the following meaings:     BondingIDs = [preID1, preID2, postID1, postID2]; Reduce = [depth1, depth2];   #
#      where depth1 sets the max depth searching from preID1 and postID1 and depth2 sets the max depth searching from preID2#
#      and postID2. NOTE ID1 and ID2 MUST be linking IDs and assumes the understanding layed out in the BondingIDs section. #
#                                                                                                                           #
#    Keep = [1, 2, 3] or Remove = [4, 5, 6]                                                                                 #
#      A string that looks like a python variable and a python list. The name Keep can be used to "keep" certain atoms when #
#      using the Reduce lists. To keep that atomIDs you desire, just add to the Keep list. Remove can "remove" certain      #
#      atomIDs, which can be used with the Reduce option or a stand alone option, if for some reason you accidentally built #
#      a reaction template with atoms where they should not belong.                                                         #
#                                                                                                                           #
#    Reduce notes:                                                                                                          #
#      The usage of the Reduce option does not check for ring opening reactions, thus it is up to the user to make sure that#
#      the depth varaible set with the Reduce list of 5 or the depth1 and depth2 variable set with the Reduce list of 2 with#
#      the coupling of BondingIDs is correct for the reaction being modeled. Setting the varaibles  write_rxn_mol2files and #
#      /or write_rxn_datafiles will allow for the visualization of the newly reduced template. In general the depth, depth1,#
#      and depth2 varaible should be large enough to make sure all bonded interations are accounted for (IE if dihedrals and#
#      impropers are being modeled the minimum value of depth or depth1 or depth2 should be 3, which finds neighbors up to  #
#      3 atoms away which is far enough for dihedral and improper interactions. Josh also recommends using the coupling of  #
#      Reduce and BondingIDs when its is able to be used due to the amount of flexibility this option allows. However, it   #
#      may not always be possible to find the "linking" atomIDs that are required for setting BondingIDs at which you can   #
#      default to using the Reduce option with list length of 5 elements. Also note that pairing any of the Reduce options  #
#      with the Keep or Remove option maybe a productive way to use the Reduce option, when dealing w/ a complex reaction.  #
#                                                                                                                           #
#    Pairing keywords and lists in a single header of the preN tagged file. BondingIDs, CreateIDs, and Reduce can be set in #
#    header as well where location doest not matter. Examples:                                                              #
#      HEADER > all2lmp: v1.4 / 9 Jan 2023 > CreateIDs = [1] BondingIDs = [7, 14, 7, 9] other text                          #
#      HEADER, BondingID = [7, 14, 7, 9] and CreateIDs = [1, 2, 3, 4] and Reduce = [3, 4]                                   #
#                                                                                                                           #
# Update generate_map_file as required.                                                                                     #
#############################################################################################################################
generate_map_file = True


#############################################################################################################################
# Options to add to molecule files (generate_map_file MUST BE True for these options to work):                              #
#   'fragment/ID/OPTION' (ordering matters and the delimiter is the '/' character).                                         #
#      index0 = 'fragment' envokes an option to add fragment section to the molecule file                                   #
#      index1 = 'ID' is user defined ID of fragment                                                                         #
#      index2 = 'OPTION' is user defined option.                                                                            #
#      Currently available options:                                                                                         #
#        'custom_charges/depth_from_edge_N' where custom_charges specifies the code to build options for bond/react's       #
#         custom_charges option and depth_from_edge_N tells the code to build a fragmentID that only has atomic charges at  #
#         a certain depth N away from any edge atom. This is useful if the charge method used could not assign a proper     #
#         charge to the edge atoms of the molecule template. This code has another option: map_near_edge_rxn_charges that   #
#         works fairly well to map charges from any dataN or postN tagged files onto edge atoms and a certain depth from    #
#         edge atoms, but this is another workaround to the charging issue of atoms near or edge atoms or edge atoms        #
#         themselves. Much like map_near_edge_rxn_charges the meaning of depth for this option is as follows:               #
#                 0 = preN tagged files will have a molecule file with a fragmentID with all atoms specified EXCEPT edge    #
#                     atoms at depth 0 from any edge.                                                                       #
#                 1 = preN tagged files will have a *.lmpmol molecule template file with a fragmentID with all atoms        #
#                     specified EXCEPT edge atoms at depth 1 from any edge.                                                 #
#                 N = preN tagged files will have a *.lmpmol molecule template file with a fragmentID with all atoms        #
#                     specified EXCEPT edge atoms at depth N from any edge.                                                 #
#         Example python string specifying this option to have a fragmentID in the preN/postN tagged molecule files that    #
#         only have atoms in the fragmentID that are 2 deep from any edge atom:                                             #
#                'fragment/example/custom_charges/depth_from_edge_2'                                                        #
#                                                                                                                           #
#   'molecule' will find moleculeIDs via a cluster analysis and order them largest-to-smallest and then adds a molecule     #
#    section to the molecule files. This info will be passed on through a list when running bond_react_merge.py in IDE      #
#    mode. Examples:                                                                                                        #
#       molecule_file_options = [] # empty python list will set no options                                                  #
#       molecule_file_options = ['molecule'] # Will add molecules section to molecule file                                  #
#       molecule_file_options = ['molecule', 'fragment/qedge/custom_charges/depth_from_edge_1']                             #
#                                                                                                                           #
# Update molecule_file_options as desired.                                                                                  #
#############################################################################################################################
molecule_file_options = ['']


#############################################################################################################################
# Option to map charges of near edge atoms from dataN or postN tagged files. *NOTE: It is assumed that there is a dataN or  #
# postN file that has an exact topological match to the preN tagged files. The dataN or postN tagged file may contain many  #
# molecules per file. Once the preN tagged file is found to be a subset of the dataN or postN tagged file the charges will  #
# be mapped up to N-neighbors away and the postN tagged files will be updated via the equivalences found via the            #
# generate_map_file option. THIS MEANS RUNNING THIS OPTION requires the generate_map_file OPTION MUST BE True. The          #
# following are valid inputs and their meanings:                                                                            #
#     False = DO NOT attempt to match topologies for preN tagged files to dataN or postN tagged files and leave charges as  #
#             is in preN/postN tagged files.                                                                                #
#     0     = match preN tagged files to dataN or postN tagged file and update only edge atom charges. This depth should be #
#             used if charge was set via the bond-increment method (int value sets depth from edge).                        #
#     1     = match preN tagged files to dataN or postN tagged file and update edge atom charges + 1 neighbor-deeper (int   #
#             value sets depth from edge).                                                                                  #
#     2     = match preN tagged files to dataN or postN tagged file and update edge atom charges + 2 neighbor-deeper (int   #
#             value sets depth from edge).                                                                                  #
#     N     = match preN tagged files to dataN or postN tagged file and update edge atom charges + N neighbor-deeper (int   #
#             value sets depth from edge).                                                                                  #
#                                                                                                                           #
# The matching seems to work well, but sometimes the code can get confused if there are symmetries on the preN tagged       #
# topology. IT IS ULTIMATELY UP TO THE USER TO CHECK IF THE CHARGES WERE MAPPED PROPERLY. The code writes comments to each  #
# atomID in the charge section of the *.moltemp file that was changed to log change for easier checking. Examples:          #
#     map_near_edge_rxn_charges = False # Will skip topology matching and charge mapping.                                   #
#     map_near_edge_rxn_charges = 0     # Will perform topology matching and charge mapping of edge atoms only              #
#     map_near_edge_rxn_charges = 1     # Will perform topology matching and charge mapping of edge atoms plus 1st neighs   #
#                                                                                                                           #
# Update map_near_edge_rxn_charges as desired (Default should be False or 0).                                               #
#############################################################################################################################
map_near_edge_rxn_charges = False


#############################################################################################################################
# Python boolean variable to assign LAMMPS's new Type Labels section(s) to the written LAMMPS data file (True or False).    #
# Examples:                                                                                                                 #
#    include_type_labels = False # Will write the "standard" format of the LAMMPS data file                                 #
#    include_type_labels = True  # Will write the Type Labels section(s) in the LAMMPS data file                            #
#                                                                                                                           #
# Update include_type_labels as desired (Default should be False).                                                          #
#############################################################################################################################
include_type_labels = True


#############################################################################################################################
# Python boolean variable to write pre/post rxn .mol2 files which ChemDraw/VMD/Avogadro can natively read (True or False).  #
# This option is useful if you would like to visualize the pre/post rxn templates in ChemDraw or VMD or Avogadro. This      #
# option also works well with generate_map_file  option since it is a good idea to check if the map file finds all          #
# equivalences properly and ChemDraw allows multiple molecules to be opened. Examples:                                      #
#     write_rxn_mol2files = True  # Will write pre/post rxn .mol2 files.                                                    #
#     write_rxn_mol2files = False # Will not write pre/post rxn .mol2 files.                                                #
#                                                                                                                           #
# Update write_rxn_mol2files as desired.                                                                                    #
#############################################################################################################################
write_rxn_mol2files = True


#############################################################################################################################
# Python boolean variable to write pre/post rxn .data files with all coeffs inside (True or False). This option is mainly   #
# for troubleshooting purposes but can be used for validation that all coeff types and all topological types are mapped     #
# properly onto the newly typed coeffs. It allows for the visualization of the *_merged.data rxn's in OVTIO, but could also #
# create files that are compatable with AutoMapper (in case users want to try generating map files using AutoMapper or use  #
# AutoMapper's reduce template functionality). Examples:                                                                    #
#     write_rxn_datafiles = True  # Will write pre/post rxn .data files.                                                    #
#     write_rxn_datafiles = False # Will not write pre/post rxn .data files.                                                #
#                                                                                                                           #
# Update write_rxn_datafiles as desired.                                                                                    #
#############################################################################################################################
write_rxn_datafiles = False


#############################################################################################################################
# Python Boolean variable (True or False) to write out a *_merged.lmpmol file for any file that has a dataN tag. This file  #
# will have all the coeff types and all topological types are mapped properly and can be used with the LAMMPS 'create_atoms'#
# command to build large random systems. Examples:                                                                          #
#     write_moleculefiles = True  # Will write dataN tagged .lmpmol files to use with LAMMPS 'create_atoms' command.        #
#     write_moleculefiles = False # Will write dataN tagged .lmpmol files to use with LAMMPS 'create_atoms' command.        #
#                                                                                                                           #
# Update write_moleculefiles as desired.                                                                                    #
#############################################################################################################################
write_moleculefiles = True


#############################################################################################################################
# Python boolean variable to print options avaible for the code at the command line interface and exit (True or False).     #
# True will print options and exit code execution and False will allow to run normally. The same printouts can be accessed  #
# with the command line override -opt or -man command line inputs. Examples                                                 #
#     print_options = True   # Will print out command line manual and exit                                                  #
#     print_options = False  # Will allow code to run to completion                                                         #
#                                                                                                                           #
# Update print_options as desired.                                                                                          #
#############################################################################################################################
print_options = False




#####################################################
### Import main from src.bond_react_merge and run ###
#####################################################
if __name__ == "__main__":  
    # Import main from src.bond_react_merge.main
    from src.bond_react_merge.main import main
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
    #       python3 bond_react_merge.py -opt                                                           #
    #              or                                                                                  #
    #       python3 bond_react_merge -man                                                              #
    #   in the terminal and the code will print out all command line option and terminate before any   #
    #   further analysis is done                                                                       #
    ####################################################################################################
    # Find command line inputs that are not the python script itself
    commandline_inputs = [i for n, i in enumerate(sys.argv) if n > 0]
    if commandline_inputs and '-gui' not in commandline_inputs: use_GUI = False
    if use_GUI or '-gui' in commandline_inputs:
        from src.bond_react_merge.GUI import bond_react_merge_GUI
        print('\n\n\nbond_react_merge is currently running in GUI mode, where all GUI inputs are intialized from bond_react_merge.\n\n\n')
        bond_react_merge_GUI(files, parent_directory, newfile, atom_style, generate_map_file, write_rxn_mol2files, write_rxn_datafiles,
                             write_moleculefiles, print_options, commandline_inputs, map_near_edge_rxn_charges, molecule_file_options,
                             include_type_labels, GUI_zoom, nfiles=maxfiles)
    else:
        # If not commandline_inputs the code is running in IDE/command line mode with out command line inputs
        # this means that all inputs to the code are handled in the Inputs section below this block of code.
        if not commandline_inputs:
            print('\n\n\nbond_react_merge is currently running in IDE mode or command line mode with no command line')
            print('inputs. This means all inputs used for the code comes from the bond_react_template_merge.py script in')
            print('the inputs section.\n\n\n')
    
        # Run main bond_react_merge function
        main(files, parent_directory, newfile, atom_style, generate_map_file, write_rxn_mol2files, write_rxn_datafiles, 
             write_moleculefiles, print_options, commandline_inputs, map_near_edge_rxn_charges, molecule_file_options, include_type_labels)