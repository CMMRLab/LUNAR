# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
September 21st, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

    *********************************************************************
    * Requirements:                                                     *
    *   python 3.7+                                                     *
    *                                                                   *
    * Dependencies:                                                     *
    *  - This code has no dependancies, but all the codes that          *
    *    this GUI can call still have their dependencies                *
    *                                                                   *
    * Run methods:                                                      *
    *   - Run from IDE                                                  *
    *   - Run from command line via:                                    *
    *       python3 LUNAR.py     # To use GUI_zoom set as the variable  *
    *       python3 LUNAR.py 125 # To set GUI_zoom at the command line  *
    *********************************************************************
"""
##############################
# Import Necessary Libraries #
##############################
from tkinter.scrolledtext import ScrolledText
from tkinter import Toplevel
import tkinter as tk
import math
import sys


#########################################################################################################
# Python int value to range from 0 to 200 to control the GUI zoom. Where it option was added to account #
# for the different ways individuals have their OS display settings setup. A GUI_zoom = 100, means that #
# default settings are used, whereas a GUI_zoom = 120 means increase GUI size by 20% and GUI_zoom = 80  #
# means decrease GUI by 20%. Examples:                                                                  #
#   GUI_zoom = 100 # use default GUI size                                                               #
#   GUI_zoom = 75  # decrease GUI size by 25%                                                           #
#                                                                                                       #
# Note that this GUI_zoom variable will be used to call all the other LUNAR codes, such that adjusting  #
# this variable will set the GUI_zoom for everything intialized from this GUI (IE atom_typing.py,       #
# all2lmp.py, bond_react_merge.py, ... etc) max have a different GUI_zoom variable then this one but    #
# when loading those GUIs from the master LUNAR GUI, this GUI_zoom variable will be used to set the GUI #
# zoom level.                                                                                           #
#                                                                                                       #
# Update GUI_zoom as desired.                                                                           #
#########################################################################################################
GUI_zoom = 100



####################
# LUNAR Master GUI #
####################
class LUNAR:
    def __init__(self, GUI_zoom):
        # Initialize window
        self.root = tk.Tk()
        self.root.title('LUNAR GUI v1.0')
        #self.root.geometry('600x400')
        self.root.resizable(width=False, height=False)
        
        # Initalize main frame
        self.frame = tk.Frame(self.root)
        self.frame.pack()
        
        # Set default sizes to use throughout this code
        try:
            from tkinter import font
            self.defaults = font.nametofont('TkTextFont').actual()
            self.font_size = self.defaults['size']
            self.font_type = self.defaults['family']
        except:
            self.font_size = 9
            self.font_type = 'Segoe UI'
            self.defaults = {'family':self.font_type, 'size':self.font_size}
        self.xpadding = 20
        self.ypadding = 10
        
        # adjust based on GUI_SF
        self.GUI_zoom = GUI_zoom
        GUI_SF = GUI_zoom/100
        font_size = int(math.ceil(GUI_SF*self.font_size))
        xpadding = int(math.ceil(GUI_SF*self.xpadding))
        ypadding = int(math.ceil(GUI_SF*self.ypadding))
        font_settings = (self.font_type, font_size)
        
        # Initalize  inputs frame
        self.inputs_frame = tk.LabelFrame(self.frame, text='Codes', font=font_settings)
        self.inputs_frame.grid(row=0, column=0, padx=xpadding, pady=ypadding)
        
        # atom_typing button
        self.atom_typing_button = tk.Button(self.inputs_frame, text='atom_typing', font=font_settings, command=self.atom_typing_GUI_in)
        self.atom_typing_button.grid(column=0, row=0)
        
        # all2lmp button
        self.all2lmp_button = tk.Button(self.inputs_frame, text='all2lmp', font=font_settings, command=self.all2lmp_GUI_in)
        self.all2lmp_button.grid(column=1, row=0)
        
        # bond_react_merge button
        self.bond_react_merge_button = tk.Button(self.inputs_frame, text='bond_react_merge', font=font_settings, command=self.bond_react_merge_GUI_in)
        self.bond_react_merge_button.grid(column=2, row=0)
        
        # bond_react_merge button
        self.bond_react_merge_prep_button = tk.Button(self.inputs_frame, text='bond_react_merge_prep', font=font_settings, command=self.bond_react_merge_prep_GUI_in)
        self.bond_react_merge_prep_button.grid(column=3, row=0)
        
        # auto_morse_bond button
        self.auto_morse_bond_button = tk.Button(self.inputs_frame, text='auto_morse_bond', font=font_settings, command=self.auto_morse_bond_GUI_in)
        self.auto_morse_bond_button.grid(column=4, row=0)
        
        # cell_builder button
        self.cell_builder_button = tk.Button(self.inputs_frame, text='cell_builder', font=font_settings, command=self.cell_builder_GUI_in)
        self.cell_builder_button.grid(column=5, row=0)
        
        # add_pi_electrons button
        self.add_pi_electrons_button = tk.Button(self.inputs_frame, text='add_pi_electrons', font=font_settings, command=self.add_pi_electrons_GUI_in)
        self.add_pi_electrons_button.grid(column=6, row=0)
        
        # auto_cluster_analysis button
        self.auto_cluster_analysis_button = tk.Button(self.inputs_frame, text='auto_cluster_analysis', font=font_settings, command=self.auto_cluster_analysis_GUI_in)
        self.auto_cluster_analysis_button.grid(column=0, row=1)
        
        # cluster_analysis button
        self.cluster_analysis_button = tk.Button(self.inputs_frame, text='cluster_analysis', font=font_settings, command=self.cluster_analysis_GUI_in)
        self.cluster_analysis_button.grid(column=1, row=1)
        
        # lmp2SYBYLmol2 button
        self.lmp2SYBYLmol2_button = tk.Button(self.inputs_frame, text='lmp2SYBYLmol2', font=font_settings, command=self.lmp2SYBYLmol2_GUI_in)
        self.lmp2SYBYLmol2_button.grid(column=2, row=1)
        
        # atom_removal button
        self.atom_removal_button = tk.Button(self.inputs_frame, text='atom_removal', font=font_settings, command=self.atom_removal_GUI_in)
        self.atom_removal_button.grid(column=3, row=1)
        
        # free_volume button
        self.free_volume_button = tk.Button(self.inputs_frame, text='free_volume', font=font_settings, command=self.free_volume_GUI_in)
        self.free_volume_button.grid(column=4, row=1)
        
        # log_analysis button
        self.log_analysis_button = tk.Button(self.inputs_frame, text='log_analysis', font=font_settings, command=self.log_analysis_GUI_in)
        self.log_analysis_button.grid(column=5, row=1)
        
        # sheet_builder button
        self.sheet_builder_button = tk.Button(self.inputs_frame, text='sheet_builder', font=font_settings, command=self.sheet_builder_GUI_in)
        self.sheet_builder_button.grid(column=6, row=1)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.inputs_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=ypadding)
        
        # Additional information
        self.about = tk.Button(self.frame, text='about', font=font_settings, command=self.about)
        self.about.grid(row=2, column=0, sticky='news', padx=xpadding, pady=int(ypadding/2))
        
        self.cite = tk.Button(self.frame, text='cite', font=font_settings, command=self.citation)
        self.cite.grid(row=3, column=0, sticky='news', padx=xpadding, pady=int(ypadding/2))
        
        #------------------------#
        # Run mainloop and close #
        #------------------------#
        self.root.protocol('WM_DELETE_WINDOW', self.closing)
        self.root.mainloop()
        
        
    #################################
    # Functions to call as commands #
    #################################
    # atom_typing GUI Run
    def atom_typing_GUI_in(self):
        # Import atom_typing variables and atom_typing_GUI
        from src.atom_typing.GUI import atom_typing_GUI
        import atom_typing as at
        
        # Run all2lmp_GUI
        atom_typing_GUI(at.topofile, at.bondfile, at.parent_directory, at.newfile, at.ff_name, at.delete_atoms, at.mass_map, at.bondorder, at.maxbonded, at.boundary,
                        at.vdw_radius_scale, at.reset_charges, at.print_options, [], at.bonds_via_distance_override, at.pdb_file, at.chargefile, at.include_comments_nta,
                        self.GUI_zoom)  
        return
    
    # all2lmp GUI Run
    def all2lmp_GUI_in(self):
        # Import all2lmp variables and all2lmp_GUI
        from src.all2lmp.GUI import all2lmp_GUI
        import all2lmp as al
        
        # Run all2lmp_GUI
        all2lmp_GUI(al.topofile, al.nta_file, al.frc_file, al.assumed, al.parent_directory, al.newfile, al.atom_style, al.ff_class,
                    al.use_auto_equivalence, al.use_assumed_auto_fill, al.reset_molids, al.reset_charges, al.write_txt_comments,
                    al.write_bond_react, al.print_options, al.use_morse_bonds, al.include_type_labels, al.add2box, al.ignore_missing_parameters,
                    al.shift, al.rotate, self.GUI_zoom, [])
        return
    
    # bond_react_merge GUI Run
    def bond_react_merge_GUI_in(self):
        # Import bond_react_merge variables and bond_react_merge_GUI
        from src.bond_react_merge.GUI import bond_react_merge_GUI
        import bond_react_merge as brm

        # Run bond_react_merge_GUI
        bond_react_merge_GUI(brm.files, brm.parent_directory, brm.newfile, brm.atom_style, brm.generate_map_file, brm.write_rxn_mol2files, 
                             brm.write_rxn_datafiles, brm.write_moleculefiles, brm.print_options, [], brm.map_near_edge_rxn_charges, brm.molecule_file_options,
                             brm.include_type_labels, self.GUI_zoom, nfiles=brm.maxfiles, scroll_bar=brm.scroll_bar)
        return
    
    # bond_react_mergeprep GUI Run
    def bond_react_merge_prep_GUI_in(self):
        # Import bond_react_merge_prep variables and bond_react_merge_prep_GUI
        from src.bond_react_merge_prep.GUI import bond_react_merge_prep_GUI
        import bond_react_merge_prep as brmp

        # Run bond_react_merge_prep_GUI
        bond_react_merge_prep_GUI(brmp.topofile, brmp.cta_file, brmp.newfile, brmp.atom_style, brmp.parent_directory, brmp.rm_unused_coeffs, self.GUI_zoom)
        return
    
    # auto_morse_bond GUI Run
    def auto_morse_bond_GUI_in(self):
        # Import auto_morse_bond variables and auto_morse_bond_GUI
        from src.auto_morse_bond.GUI import auto_morse_bond_GUI
        #from src.auto_morse_bond.GUI_simple import auto_morse_bond_GUI
        import auto_morse_bond_update as amb

        # Run auto_morse_bond_GUI
        auto_morse_bond_GUI(amb.topofile, amb.morsefile, amb.parent_directory, amb.newfile, amb.mass_map, amb.min_bond_length, amb.coeffs2skip,
                            amb.radius_specs, amb.alpha_specs, amb.alpha_scale, amb.files2write, amb.atom_style, amb.zero_effected_xterms,
                            amb.bondbreak_scale, amb.ff_class, amb.include_type_labels, amb.class2xe_update, amb.include_rcut, self.GUI_zoom)
        return
    
    # cell_builder GUI Run
    def cell_builder_GUI_in(self):
        # Import cell_builder variables and cell_builder_GUI
        from src.cell_builder.GUI import cell_builder_GUI
        import cell_builder as cb

        # Run cell_builder_GUI
        cell_builder_GUI(cb.files, cb.force_field_joining, cb.duplicate, cb.distance_scale, cb.newfile, cb.atom_style, cb.parent_directory, cb.max_rotations, cb.reset_molids,
                         cb.unwrap_atoms_via_image_flags, cb.include_type_labels, cb.group_monomers_locally, cb.seed, cb.domain, cb.maxtry, cb.tolerance, cb.mixing_rule,
                         cb.boundary, self.GUI_zoom, nfiles=cb.maxfiles, scroll_bar=cb.scroll_bar)
        return
    
    # add_pi_electrons GUI Run
    def add_pi_electrons_GUI_in(self):
        # Import add_pi_electrons variables and add_pi_electrons_GUI
        from src.add_pi_electrons.GUI import add_pi_electrons_GUI 
        import add_pi_electrons as a

        # Run add_pi_electrons_GUI
        add_pi_electrons_GUI(a.topofile, a.types2convert, a.atom_style, a.reset_charges, a.net_zero_charge, a.convert2cg1, a.add_pi_electrons,
                             a.parent_directory, a.newfile, a.include_type_labels, a.neighbor_charge_constraint, a.reset_simulation_cell, self.GUI_zoom)
        return
    
    # auto_cluster_analysis GUI Run
    def auto_cluster_analysis_GUI_in(self):
        # Import auto_cluster_analysis
        from src.auto_clusters_GUI import GUI
        import auto_cluster_analysis as ac

        # Run auto_clusters_GUI
        GUI(ac.files_directory, ac.N0, ac.fav, ac.txtfile, ac.newfile, self.GUI_zoom)
        return
    
    # cluster_analysis GUI Run
    def cluster_analysis_GUI_in(self):
        # Import cluster_analysis
        from src.clusters_GUI import GUI
        import cluster_analysis as ca

        # Run clusters_GUI
        GUI(ca.topofile, ca.N0, ca.txtfile, ca.fav, self.GUI_zoom, pflag=True)
        return
    
    # lmp2SYBYLmol2 GUI Run
    def lmp2SYBYLmol2_GUI_in(self):
        # Import lmp2SYBYLmol2
        from src.lmp2SYBYLmol2_GUI import GUI
        import lmp2SYBYLmol2 as l

        # Run lmp2SYBYLmol2_GUI
        GUI(l.topofile, l.parent_directory, l.remove_PBC_bonds, l.mass_map, self.GUI_zoom)
        return
    
    # atom_removal GUI Run
    def atom_removal_GUI_in(self):
        # Import atom_removal
        from src.atom_removal.GUI import atom_removal_GUI
        import atom_removal as ar

        # Run atom_removal_GUI
        atom_removal_GUI(ar.topofile, ar.parent_directory, ar.newfile, ar.atom_style, ar.atoms2remove, ar.include_type_labels, self.GUI_zoom, method=ar.method)
        return
    
    # free_volume GUI Run
    def free_volume_GUI_in(self):
        # Import free_volume
        from src.free_volume.GUI import free_volume_GUI
        import free_volume as fv

        # Run free_volume GUI
        free_volume_GUI(fv.topofile, fv.max_voxel_size, fv.mass_map, fv.vdw_radius, fv.boundary, fv.parent_directory,
                        fv.compute_free_volume_distributions, fv.files2write, fv.run_mode, fv.probe_diameter, fv.vdw_method,
                        fv.CUDA_threads_per_block_atoms, fv.CUDA_threads_per_block_voxels, self.GUI_zoom)
        return
    
    # log_analysis GUI Run
    def log_analysis_GUI_in(self):
        # Import log_analysis
        from src.log_analysis.GUI import GUI
        import log_analysis as l

        # Run GUI
        GUI(l.settings, self.GUI_zoom)
        return
    
    # sheet_builder GUI Run
    def sheet_builder_GUI_in(self):
        # Import sheet_builder
        from src.sheet_builder.GUI import sheet_builder_GUI
        import sheet_builder as sb

        # Run GUI
        sheet_builder_GUI(sb.sheet_basename, sb.symmetric_tube_basename, sb.chiral_tube_basename, sb.run_mode, sb.parent_directory, sb.length_in_perpendicular, sb.length_in_edgetype,
                          sb.sheet_edgetype, sb.types, sb.bond_length, sb.sheet_layer_spacing, sb.sheet_nlayers, sb.stacking, sb.plane, sb.tube_edgetype, sb.tube_layer_spacing,
                          sb.symmetric_ntubes, sb.symmetric_length, sb.diameter, sb.n, sb.m, sb.chiral_length, sb.symmetric_tube_axis, sb.chiral_tube_axis, sb.find_bonds, sb.periodic_bonds,
                          sb.charges, sb.masses, sb.functional_seed, sb.functional_atoms, sb.terminating_atoms, self.GUI_zoom)
        return
    
    # Function to pop-up scrollable text
    def popup(self, out, title='Outputs', width=150, height=30):
        page = Toplevel(self.root)
        page.title(title)
        outputs = ScrolledText(page, height=height, width=width, font=('consolas', '12', 'normal'))
        outputs.pack()
        outputs.insert(tk.INSERT, '\n'.join(out))
        outputs.config(state=tk.DISABLED)
        page.resizable(width=False, height=False)
        return
    
    
    # Def citations page
    def citation(self):
        txt = ['']
        txt.append('LUNAR: Automated Input Generation and Analysis for Reactive LAMMPS Simulations')
        txt.append('at https://doi.org/10.1021/acs.jcim.4c00730')
        self.popup(txt, title='cite', width=80, height=5)
        return
        
    # Def about page
    def about(self):
        txt = []
        txt.append('LUNAR is a set of codes developed by Josh Kemppainen during his Ph.D. at Michigan Technological')
        txt.append('University. The LUNAR codes aim at simplifying the MD model-building process required to use')
        txt.append('LAMMPS and features found within LAMMPS such as fix bond/react (aka REACTER). The LUNAR codes')
        txt.append('allow for a heuristic method of learning the basics of MD inputs and offer the ability to change')
        txt.append('force fields from bond order -> fix bond or fix bond -> fix bond or fix bond ->  bond order')
        txt.append('within a single consistent set of codes.')
        txt.append('')
        txt.append('LUNAR is a medium-difficulty method of building inputs to LAMMPS which makes it suitable for') 
        txt.append('beginners learning MD and LAMMPS, but also offers the flexibility for more experienced')
        txt.append('individuals to build models exactly as they wish.')
        txt.append('')
        txt.append('LUNAR is coded entirely in Python v3 and all codes assume Python 3.7 or greater. If you need')
        txt.append('assistance learning LUNAR or would like more features added to LUNAR reach out to Josh at:')
        txt.append('    jdkemppa@mtu.edu')
        txt.append('')
        txt.append('If you are a developer and would like to add a package to LUNAR you may wish to do so yourself')
        txt.append('or get in contact with Josh to discuss the best integration method. It would be great to add')
        txt.append('more charging methods built into LUNAR and also allow for a wider range of force fields, so if')
        txt.append('If you have expertise in these fields and would like a collaboration please reach out to Josh at:')
        txt.append('    jdkemppa@mtu.edu')
        self.popup(txt, title='about', width=100, height=25)


    
    # Closing command    
    def closing(self):
        print('Terminating LUNAR GUI'); self.root.destroy();
        # if tk.messagebox.askyesno(title='Quit?', message='Do you really want to quit?'):
        #     print('Terminating LUNAR GUI'); self.root.destroy();
        return
        
        
########################
# Run LUNAR Master GUI #
########################
if __name__ == "__main__": 
    commandline_inputs = [i for n, i in enumerate(sys.argv) if n > 0]
    if len(commandline_inputs) == 1: 
        GUI_zoom = int(commandline_inputs[0])
        print('Updating GUI_zoom from command line as: ', GUI_zoom)
    LUNAR(GUI_zoom)