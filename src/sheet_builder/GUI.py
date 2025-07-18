# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
September 21st, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.GUI_scale_settings as GUI_scale_settings
from src.sheet_builder.main import main
import src.io_functions as io_functions
import src.py_script_modifier as psm
from tkinter.scrolledtext import ScrolledText
from tkinter import filedialog
from tkinter import Toplevel
from tkinter import ttk
import tkinter as tk
import threading
import traceback
import math
import os


##########################
# LUNAR/sheet_builder GUI #
##########################
class sheet_builder_GUI:
    def __init__(self, sheet_basename, symmetric_tube_basename, chiral_tube_basename, run_mode, parent_directory, length_in_perpendicular, length_in_edgetype,
                 sheet_edgetype, types, bond_length, sheet_layer_spacing, sheet_nlayers, stacking, plane, tube_edgetype, tube_layer_spacing,
                 symmetric_ntubes, symmetric_length, diameter, n, m, chiral_length, symmetric_tube_axis, chiral_tube_axis, find_bonds, periodic_bonds,
                 charges, masses, seed, functional_atoms, terminating_atoms, grafting_files, minimum_distance, cutter, GUI_zoom):
        
        # Find present working directory
        self.pwd = os.getcwd()
        self.filename = os.path.join(self.pwd, 'sheet_builder.py')
        self.charges = charges
        self.masses = masses

        # Initialize window
        self.root = tk.Tk()
        self.root.title('LUNAR/sheet_builder.py GUI v1.0')
        self.root.resizable(width=False, height=False)
        
        
        # Setup tabs
        self.tabControl = ttk.Notebook(self.root)
        self.tab1 = ttk.Frame(self.tabControl)
        self.tab2 = ttk.Frame(self.tabControl)
        self.tabControl.add(self.tab1, text='Main and Sheets mode')
        self.tabControl.add(self.tab2, text='Symmetric and Chiral tube modes')
        self.tabControl.pack(expand=1, fill='both')
        #self.root.geometry('600x400')
        
        # Initalize main frame
        self.frame1 = tk.Frame(self.tab1)
        self.frame1.pack()
        self.frame2 = tk.Frame(self.tab2)
        self.frame2.pack()
        
        #-----------------------------------------------#
        # Set default sizes to use throughout this code #
        #-----------------------------------------------#
        # Set/get defaults
        try:
            from tkinter import font
            self.defaults = font.nametofont('TkTextFont').actual()
            self.font_size = self.defaults['size']
            self.font_type = self.defaults['family']
        except:
            self.font_size = 9
            self.font_type = 'Segoe UI'
            self.defaults = {'family':self.font_type, 'size':self.font_size}
            
        # Check if user specified any other font settings
        font_settings = GUI_scale_settings.font_settings
        self.int_font = font.nametofont("TkDefaultFont") # font parameters used in internal TK dialogues
        self.icon_font = font.nametofont("TkIconFont") # font used in file selection dialogue
        if 'size' in font_settings:
            if isinstance(font_settings['size'], (int, float)):
               self.font_size = font_settings['size'] 
        if 'type' in font_settings:
            if isinstance(font_settings['type'], str):
               self.font_type = font_settings['type'] 
        if 'dialog_size' in font_settings:
            if isinstance(font_settings['dialog_size'], (int, float)):
               self.int_font.configure(size=font_settings['dialog_size'])
               self.icon_font.configure(size=font_settings['dialog_size'])
        if 'dialog_type' in font_settings:
            if isinstance(font_settings['dialog_type'], str):
               self.int_font.configure(family=font_settings['dialog_type'])
               self.icon_font.configure(family=font_settings['dialog_type'])
        self.defaults = {'family':self.font_type, 'size':self.font_size}

        self.xpadding = 20
        self.ypadding = 10
        self.maxwidth = 130

        # Check if user specified any other nong-global scaling settings
        scale_settings = GUI_scale_settings.screen_settings
        if 'scaling_factor' in scale_settings:
            if isinstance(scale_settings['scaling_factor'], (int, float)):
               self.xpadding = int(self.xpadding/scale_settings['scaling_factor'])
               self.ypadding = int(self.ypadding/scale_settings['scaling_factor'])
               self.maxwidth = int(self.maxwidth/scale_settings['scaling_factor'])
        
        # adjust based on GUI_SF
        GUI_SF = GUI_zoom/100
        font_size = int(math.ceil(GUI_SF*self.font_size))
        xpadding = int(math.ceil(GUI_SF*self.xpadding))
        ypadding = int(math.ceil(GUI_SF*self.ypadding))
        maxwidth = int(math.ceil(GUI_SF*self.maxwidth))
        self.font_settings = (self.font_type, font_size)
        self.GUI_zoom = GUI_zoom
        
        # Set up types to load in typeN drop down menu. If a type is provide by the user that is not in defaults, append user type
        self.types = ['', 'C', 'B', 'N', 'cp', 'cg1', 'cg1|cge', 'bbn', 'bbn|bbe', 'nbn', 'nbn|nbe']
        for i in types:
            if types[i] not in self.types: 
                self.types.append(types[i])
                
        #--------------#
        # Inputs frame #
        #--------------#
        # Initalize  inputs frame
        self.inputs_frame = tk.LabelFrame(self.frame1, text='Global Inputs', font=self.font_settings)
        self.inputs_frame.grid(row=0, column=0, columnspan=2, padx=xpadding, pady=ypadding)

        # parent_directory entry
        self.parent_directory = tk.Entry(self.inputs_frame, width=int(1.0*maxwidth), font=self.font_settings)
        self.parent_directory.insert(0, parent_directory)
        self.parent_directory.grid(column=1, row=0, columnspan=8)
        self.dir_button = tk.Button(self.inputs_frame, text='parent_directory', font=self.font_settings, command=self.directory_path)
        self.dir_button.grid(column=0, row=0)
        
        # type1 drop down menu
        self.type1 = ttk.Combobox(self.inputs_frame, values=self.types, width=int(maxwidth/12.75), font=self.font_settings)
        self.type1.current(self.types.index(types[1]))
        self.type1.grid(column=1, row=1)
        self.type1_label = tk.Label(self.inputs_frame, text='type1', font=self.font_settings)
        self.type1_label.grid(column=0, row=1)
        
        # type2 entry
        self.type2 = ttk.Combobox(self.inputs_frame, values=self.types, width=int(maxwidth/12.75), font=self.font_settings)
        self.type2.current(self.types.index(types[2]))
        self.type2.grid(column=3, row=1)
        self.type2_label = tk.Label(self.inputs_frame, text='type2', font=self.font_settings)
        self.type2_label.grid(column=2, row=1)
        
        # type3 entry
        self.type3 = ttk.Combobox(self.inputs_frame, values=self.types, width=int(maxwidth/12.75), font=self.font_settings)
        self.type3.current(self.types.index(types[3]))
        self.type3.grid(column=5, row=1)
        self.type3_label = tk.Label(self.inputs_frame, text='type3', font=self.font_settings)
        self.type3_label.grid(column=4, row=1)
        
        # type4 entry
        self.type4 = ttk.Combobox(self.inputs_frame, values=self.types, width=int(maxwidth/12.75), font=self.font_settings)
        self.type4.current(self.types.index(types[4]))
        self.type4.grid(column=7, row=1)
        self.type4_label = tk.Label(self.inputs_frame, text='type4', font=self.font_settings)
        self.type4_label.grid(column=6, row=1)
        
        # bond_length entry
        self.bond_length = tk.Entry(self.inputs_frame, width=int(maxwidth/10), font=self.font_settings)
        self.bond_length.insert(0, bond_length)
        self.bond_length.grid(column=1, row=2)
        self.bond_length_label = tk.Label(self.inputs_frame, text='bond_length (Å)', font=self.font_settings)
        self.bond_length_label.grid(column=0, row=2)
        
        # find_bonds drop down menu
        styles = [True, False]
        self.find_bonds = ttk.Combobox(self.inputs_frame, values=styles, width=int(maxwidth/12.75), font=self.font_settings)
        self.find_bonds.current(styles.index(find_bonds))
        self.find_bonds.grid(column=3, row=2)
        self.find_bonds_label = tk.Label(self.inputs_frame, text='find_bonds', font=self.font_settings)
        self.find_bonds_label.grid(column=2, row=2)
        
        # periodic_bonds drop down menu
        styles = [True, False]
        self.periodic_bonds = ttk.Combobox(self.inputs_frame, values=styles, width=int(maxwidth/12.75), font=self.font_settings)
        self.periodic_bonds.current(styles.index(periodic_bonds))
        self.periodic_bonds.grid(column=5, row=2)
        self.periodic_bonds_label = tk.Label(self.inputs_frame, text='periodic_bonds', font=self.font_settings)
        self.periodic_bonds_label.grid(column=4, row=2)
        
        # functional seed entry
        self.seed = tk.Entry(self.inputs_frame, width=int(maxwidth/10), font=self.font_settings)
        self.seed.insert(0, seed)
        self.seed.grid(column=7, row=2)
        self.seed_label = tk.Label(self.inputs_frame, text='seed', font=self.font_settings)
        self.seed_label.grid(column=6, row=2)
        
        # terminating_atoms entry
        self.terminating_atoms = tk.Entry(self.inputs_frame, width=int(0.55*maxwidth), font=self.font_settings)
        self.terminating_atoms.insert(0, terminating_atoms)
        self.terminating_atoms.grid(column=1, row=3, columnspan=4)
        self.terminating_atoms_label = tk.Label(self.inputs_frame, text='terminating_atoms', font=self.font_settings)
        self.terminating_atoms_label.grid(column=0, row=3)
        
        # minimum_distance entry
        self.minimum_distance = tk.Entry(self.inputs_frame, width=int(0.25*maxwidth), font=self.font_settings)
        self.minimum_distance.insert(0, minimum_distance)
        self.minimum_distance.grid(column=6, row=3, columnspan=2)
        self.minimum_distance_label = tk.Label(self.inputs_frame, text='minimum_distance', font=self.font_settings)
        self.minimum_distance_label.grid(column=5, row=3)
        
        
        # functional_atoms entry
        self.functional_atoms = tk.Entry(self.inputs_frame, width=int(0.99*maxwidth), font=self.font_settings)
        self.functional_atoms.insert(0, functional_atoms)
        self.functional_atoms.grid(column=1, row=4, columnspan=7)
        self.functional_atoms_label = tk.Label(self.inputs_frame, text='functional_atoms', font=self.font_settings)
        self.functional_atoms_label.grid(column=0, row=4)
        
        # grafting_files entry
        self.grafting_files = tk.Entry(self.inputs_frame, width=int(0.99*maxwidth), font=self.font_settings)
        self.grafting_files.insert(0, grafting_files)
        self.grafting_files.grid(column=1, row=5, columnspan=7)
        self.grafting_files_label = tk.Label(self.inputs_frame, text='grafting_files', font=self.font_settings)
        self.grafting_files_label.grid(column=0, row=5)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.inputs_frame.winfo_children():
            widget.grid_configure(padx=int(xpadding/2), pady=int(ypadding/1))
            
        
        #--------------#
        # Cutter frame #
        #--------------#
        # Initalize cutter frame
        self.cutter_frame = tk.LabelFrame(self.frame1, text='Cutter', font=self.font_settings)
        self.cutter_frame.grid(row=1, column=0, columnspan=2, sticky='news', padx=self.xpadding, pady=self.ypadding)
        
        self.cutter = cutter
        self.ncuts = len(self.cutter)
        self.cutter = []
        for n in range(1, self.ncuts+1):
            cut = cutter[n-1]
            
            self.cut = tk.Entry(self.cutter_frame, width=int(1.125*self.maxwidth), font=self.font_settings)
            self.cut.grid(column=0, row=n, columnspan=3)
            self.cut.insert(0, cut)
            self.cutter.append(self.cut)
            
        # Button to add a file
        self.add_button = tk.Button(self.cutter_frame, text='add cutter to stack', font=self.font_settings, width=int(self.maxwidth/3), command=self.add2stack)
        self.add_button.grid(column=0, row=self.ncuts+1, columnspan=1)
            
        # Button to remove a file
        self.remove_button = tk.Button(self.cutter_frame, text='remove last cutter from stack', font=self.font_settings, width=int(self.maxwidth/3.5), command=self.remove_last)
        self.remove_button.grid(column=1, row=self.ncuts+1, sticky='news', columnspan=1)
        
        # Button to clear all files
        self.clear_button = tk.Button(self.cutter_frame, text='clear stack', font=self.font_settings, width=int(self.maxwidth/3.5), command=self.clear_all)
        self.clear_button.grid(column=2, row=self.ncuts+1, columnspan=1)
            
            
        # Add padding to all frames in self.cutter_frame
        for widget in self.cutter_frame.winfo_children():
            widget.grid_configure(padx=int(xpadding/2), pady=int(ypadding/1))
            
        
        #--------------#
        # Sheets frame #
        #--------------#
        # Initalize sheets frame
        self.sheets_frame = tk.LabelFrame(self.frame1, text='Sheets mode', font=self.font_settings)
        self.sheets_frame.grid(row=2, column=0, columnspan=2, sticky='news', padx=xpadding, pady=int(maxwidth/9.5))
        
        # sheet_basename entry
        self.sheet_basename = tk.Entry(self.sheets_frame, width=maxwidth, font=self.font_settings)
        self.sheet_basename.insert(0, sheet_basename)
        self.sheet_basename.grid(column=1, row=0, columnspan=8)
        self.sheet_basename_label = tk.Label(self.sheets_frame, text='basename', font=self.font_settings)
        self.sheet_basename_label.grid(column=0, row=0)
        
        # length_in_edgetype entry
        self.length_in_edgetype = tk.Entry(self.sheets_frame, width=int(maxwidth/10), font=self.font_settings)
        self.length_in_edgetype.insert(0, length_in_edgetype)
        self.length_in_edgetype.grid(column=1, row=1)
        self.length_in_edgetype_label = tk.Label(self.sheets_frame, text='edge length (Å)', font=self.font_settings)
        self.length_in_edgetype_label.grid(column=0, row=1)
        
        # length_in_perpendicular entry
        self.length_in_perpendicular = tk.Entry(self.sheets_frame, width=int(maxwidth/10), font=self.font_settings)
        self.length_in_perpendicular.insert(0, length_in_perpendicular)
        self.length_in_perpendicular.grid(column=3, row=1)
        self.length_in_perpendicular_label = tk.Label(self.sheets_frame, text='perpendicular length (Å)', font=self.font_settings)
        self.length_in_perpendicular_label.grid(column=2, row=1)
        
        # sheet_layer_spacing entry
        self.sheet_layer_spacing = tk.Entry(self.sheets_frame, width=int(maxwidth/10), font=self.font_settings)
        self.sheet_layer_spacing.insert(0, sheet_layer_spacing)
        self.sheet_layer_spacing.grid(column=5, row=1)
        self.sheet_layer_spacing_label = tk.Label(self.sheets_frame, text='layer spacing (Å)', font=self.font_settings)
        self.sheet_layer_spacing_label.grid(column=4, row=1)
        
        # sheet_nlayers entry
        self.sheet_nlayers = tk.Entry(self.sheets_frame, width=int(maxwidth/12), font=self.font_settings)
        self.sheet_nlayers.insert(0, sheet_nlayers)
        self.sheet_nlayers.grid(column=7, row=1)
        self.sheet_nlayers_label = tk.Label(self.sheets_frame, text='nlayers', font=self.font_settings)
        self.sheet_nlayers_label.grid(column=6, row=1)

        # sheet_edgetype drop down
        styles = ['armchair', 'zigzag']
        self.sheet_edgetype = ttk.Combobox(self.sheets_frame, values=styles, width=int(maxwidth/12.5), font=self.font_settings)
        self.sheet_edgetype.current(styles.index(sheet_edgetype))
        self.sheet_edgetype.grid(column=1, row=2)
        self.sheet_edgetype_label = tk.Label(self.sheets_frame, text='edge type', font=self.font_settings)
        self.sheet_edgetype_label.grid(column=0, row=2)
        
        # stacking drop down
        styles = ['AA', 'AB', 'ABC']
        self.stacking = ttk.Combobox(self.sheets_frame, values=styles, width=int(maxwidth/12.5), font=self.font_settings)
        self.stacking.current(styles.index(stacking))
        self.stacking.grid(column=3, row=2)
        self.stacking_label = tk.Label(self.sheets_frame, text='stacking', font=self.font_settings)
        self.stacking_label.grid(column=2, row=2)
        
        # stacking drop down
        styles = ['xy', 'xz', 'yz']
        self.plane = ttk.Combobox(self.sheets_frame, values=styles, width=int(maxwidth/12.5), font=self.font_settings)
        self.plane.current(styles.index(plane))
        self.plane.grid(column=5, row=2)
        self.plane_label = tk.Label(self.sheets_frame, text='plane', font=self.font_settings)
        self.plane_label.grid(column=4, row=2)
        
        # Button to run in 'sheet' mode
        self.run = tk.Button(self.sheets_frame, text="Run LUNAR/sheet_builder.py in 'sheet' mode", font=self.font_settings, command=self.run_in_sheet_mode)
        self.run.grid(row=3, column=0, columnspan=8, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.sheets_frame.winfo_children():
            widget.grid_configure(padx=int(xpadding/2), pady=int(ypadding/1))
            
            
        #----------------------#
        # symmetric tube frame #
        #----------------------#
        # Initalize symmetric tubes frame
        self.symmetric_tube_frame = tk.LabelFrame(self.frame2, text='Symmetric tubes mode', font=self.font_settings)
        self.symmetric_tube_frame.grid(row=2, column=0, columnspan=2, sticky='news', padx=xpadding, pady=int(maxwidth/9.5))
        
        # symmetric_tube_basename entry
        self.symmetric_tube_basename = tk.Entry(self.symmetric_tube_frame, width=maxwidth, font=self.font_settings)
        self.symmetric_tube_basename.insert(0, symmetric_tube_basename)
        self.symmetric_tube_basename.grid(column=1, row=0, columnspan=8)
        self.symmetric_tube_basename_label = tk.Label(self.symmetric_tube_frame, text='basename', font=self.font_settings)
        self.symmetric_tube_basename_label.grid(column=0, row=0)
        
        # symmetric_length entry
        self.symmetric_length = tk.Entry(self.symmetric_tube_frame, width=int(maxwidth/10), font=self.font_settings)
        self.symmetric_length.insert(0, symmetric_length)
        self.symmetric_length.grid(column=1, row=1)
        self.symmetric_length_label = tk.Label(self.symmetric_tube_frame, text='tube length (Å)', font=self.font_settings)
        self.symmetric_length_label.grid(column=0, row=1)
        
        # diameter entry
        self.diameter = tk.Entry(self.symmetric_tube_frame, width=int(maxwidth/10), font=self.font_settings)
        self.diameter.insert(0, diameter)
        self.diameter.grid(column=3, row=1)
        self.diameter_label = tk.Label(self.symmetric_tube_frame, text='tube diameter (Å)', font=self.font_settings)
        self.diameter_label.grid(column=2, row=1)
        
        # tube_layer_spacing entry
        self.tube_layer_spacing = tk.Entry(self.symmetric_tube_frame, width=int(maxwidth/10), font=self.font_settings)
        self.tube_layer_spacing.insert(0, tube_layer_spacing)
        self.tube_layer_spacing.grid(column=5, row=1)
        self.tube_layer_spacing_label = tk.Label(self.symmetric_tube_frame, text='layer spacing (Å)', font=self.font_settings)
        self.tube_layer_spacing_label.grid(column=4, row=1)
        
        # symmetric_ntubes entry
        self.symmetric_ntubes = tk.Entry(self.symmetric_tube_frame, width=int(maxwidth/10), font=self.font_settings)
        self.symmetric_ntubes.insert(0, symmetric_ntubes)
        self.symmetric_ntubes.grid(column=7, row=1)
        self.symmetric_ntubes_label = tk.Label(self.symmetric_tube_frame, text='ntubes', font=self.font_settings)
        self.symmetric_ntubes_label.grid(column=6, row=1)
        
        # tube_edgetype drop down
        styles = ['armchair', 'zigzag']
        self.tube_edgetype = ttk.Combobox(self.symmetric_tube_frame, values=styles, width=int(maxwidth/12.5), font=self.font_settings)
        self.tube_edgetype.current(styles.index(tube_edgetype))
        self.tube_edgetype.grid(column=1, row=2)
        self.tube_edgetype_label = tk.Label(self.symmetric_tube_frame, text='edge type', font=self.font_settings)
        self.tube_edgetype_label.grid(column=0, row=2)
        
        # symmetric_tube_axis drop down
        styles = ['x', 'y', 'z']
        self.symmetric_tube_axis = ttk.Combobox(self.symmetric_tube_frame, values=styles, width=int(maxwidth/12.5), font=self.font_settings)
        self.symmetric_tube_axis.current(styles.index(symmetric_tube_axis))
        self.symmetric_tube_axis.grid(column=3, row=2)
        self.symmetric_tube_axis_label = tk.Label(self.symmetric_tube_frame, text='axis', font=self.font_settings)
        self.symmetric_tube_axis_label.grid(column=2, row=2)
        
        # Button to run in 'symmetric-tube' mode
        self.run = tk.Button(self.symmetric_tube_frame, text="Run LUNAR/sheet_builder.py in 'symmetric-tube' mode", font=self.font_settings, command=self.run_in_symmetric_tube_mode)
        self.run.grid(row=3, column=0, columnspan=8, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.symmetric_tube_frame.winfo_children():
            widget.grid_configure(padx=int(xpadding/2), pady=int(ypadding/1))
            
        #-------------------#
        # Chiral tube frame #
        #-------------------#
        # Initalize chiral tubes frame
        self.chiral_tube_frame = tk.LabelFrame(self.frame2, text='Chiral tubes mode', font=self.font_settings)
        self.chiral_tube_frame.grid(row=3, column=0, columnspan=2, sticky='news', padx=xpadding, pady=int(maxwidth/9.5))
        
        # symmetric_tube_basename entry
        self.chiral_tube_basename = tk.Entry(self.chiral_tube_frame, width=maxwidth, font=self.font_settings)
        self.chiral_tube_basename.insert(0, chiral_tube_basename)
        self.chiral_tube_basename.grid(column=1, row=0, columnspan=8)
        self.chiral_tube_basename_label = tk.Label(self.chiral_tube_frame, text='basename', font=self.font_settings)
        self.chiral_tube_basename_label.grid(column=0, row=0)
        
        # chiral_length entry
        self.chiral_length = tk.Entry(self.chiral_tube_frame, width=int(maxwidth/10), font=self.font_settings)
        self.chiral_length.insert(0, chiral_length)
        self.chiral_length.grid(column=1, row=1)
        self.chiral_length_label = tk.Label(self.chiral_tube_frame, text='tube length (Å)', font=self.font_settings)
        self.chiral_length_label.grid(column=0, row=1)
        
        # n entry
        self.n = tk.Entry(self.chiral_tube_frame, width=int(maxwidth/10), font=self.font_settings)
        self.n.insert(0, n)
        self.n.grid(column=3, row=1)
        self.n_label = tk.Label(self.chiral_tube_frame, text='n', font=self.font_settings)
        self.n_label.grid(column=2, row=1)
        
        # m entry
        self.m = tk.Entry(self.chiral_tube_frame, width=int(maxwidth/10), font=self.font_settings)
        self.m.insert(0, m)
        self.m.grid(column=5, row=1)
        self.m_label = tk.Label(self.chiral_tube_frame, text='m', font=self.font_settings)
        self.m_label.grid(column=4, row=1)
        
        # symmetric_tube_axis drop down
        styles = ['x', 'y', 'z']
        self.chiral_tube_axis = ttk.Combobox(self.chiral_tube_frame, values=styles, width=int(maxwidth/12.5), font=self.font_settings)
        self.chiral_tube_axis.current(styles.index(chiral_tube_axis))
        self.chiral_tube_axis.grid(column=7, row=1)
        self.chiral_tube_axis_label = tk.Label(self.chiral_tube_frame, text='axis', font=self.font_settings)
        self.chiral_tube_axis_label.grid(column=6, row=1)
        
        # Button to run in 'chiral-tube' mode
        self.run = tk.Button(self.chiral_tube_frame, text="Run LUNAR/sheet_builder.py in 'chiral-tube' mode", font=self.font_settings, command=self.run_in_chiral_tube_mode)
        self.run.grid(row=2, column=0, columnspan=8, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.chiral_tube_frame.winfo_children():
            widget.grid_configure(padx=int(xpadding/2), pady=int(ypadding/1))

    
        #-----------------#
        # update defaults #
        #-----------------#
        self.update = tk.Button(self.frame1, text='Save the current GUI settings as the default GUI settings', font=self.font_settings, command=self.update_py_script)
        self.update.grid(row=3, column=0, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
        self.update = tk.Button(self.frame2, text='Save the current GUI settings as the default GUI settings', font=self.font_settings, command=self.update_py_script)
        self.update.grid(row=5, column=0, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
        #------------#
        # Quick help #
        #------------#
        self.quick_help = tk.Button(self.frame1, text='Quick help', font=self.font_settings, command=self.quickhelp)
        self.quick_help.grid(row=3, column=1, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
        self.quick_help = tk.Button(self.frame2, text='Quick help', font=self.font_settings, command=self.quickhelp)
        self.quick_help.grid(row=5, column=1, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
        #------------------------#
        # Run mainloop and close #
        #------------------------#
        self.root.protocol('WM_DELETE_WINDOW', self.closing)
        self.root.mainloop()
        
        
        
    #################################
    # Functions to call as commands #
    #################################
    # Quick help button
    def quickhelp(self):
        try: # Try to get text from GUI_help_page.txt file
            txt = os.path.join(self.pwd, 'src/GUI_quick_help_pages/sheet_builder.txt')
            logged = []
            with open(txt, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    # Strip comment's and split by whitespace
                    line = line.split('#')[0]
                    line = line.rstrip()
                    logged.append(line)
        except: # except something failed
            logged.append('FAILED to read LUNAR/src/GUI_quick_help_pages/sheet_builder.txt document.')
            logged.append('Most likely cause is the sheet_builder.txt file was renamed in the')
            logged.append('LUNAR/src/GUI_quick_help_pages/ directory or directory names were changed.')
        self.popup(logged, title='Quick help')
        return   
    
    # Function to add to analysis options
    def add2stack(self):
        try: self.add_overloaded_analysis()
        except: print('GUI failed to add additionaly analysis to stack')
        return
    
    # Function to add files to GUI, during overload conditions
    def add_overloaded_analysis(self):    
        # adjust based on GUI_SF
        GUI_SF = self.GUI_zoom/100
        xpadding = int(math.ceil(GUI_SF*self.xpadding))
        ypadding = int(math.ceil(GUI_SF*self.ypadding))
        
        # Add file box
        self.ncuts += 1        
        self.cut = tk.Entry(self.cutter_frame, width=int(1.125*self.maxwidth), font=self.font_settings)
        self.cut.grid(column=0, row=self.ncuts, columnspan=3)
        self.cut.insert(0, '')
        self.cutter.append(self.cut)
        
        # adjust packing of other things in inputs frame
        self.add_button.grid(column=0, row=self.ncuts+1, columnspan=1)
        self.remove_button.grid(column=1, row=self.ncuts+1, sticky='news', columnspan=1)
        self.clear_button.grid(column=2, row=self.ncuts+1, columnspan=1)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.cutter_frame.winfo_children():
            widget.grid_configure(padx=int(xpadding/2), pady=int(ypadding/1))
        return
    
    # Function to remove last anaylsis
    def remove_last(self):
        used = []
        for n, i in enumerate(self.cutter):
            cut = i.get()
            if cut != '': used.append(n)
        if used:
            last_add = max(used)
            self.cutter[last_add].delete(0, tk.END)
        else: print('No files or tags left to remove')
        return
    
    # Function to clear all files
    def clear_all(self):
        for n, i in enumerate(self.cutter):
            try: self.cutter[n].delete(0, tk.END)
            except: pass
        return
    
    # Function to get directory
    def directory_path(self):
        path = filedialog.askdirectory()
        if path:
            path = os.path.relpath(path)
            self.parent_directory.delete(0, tk.END); self.parent_directory.insert(0, path);
        return
    
    # Closing command    
    def closing(self):
        print('Terminating sheet removal GUI'); self.root.destroy();
        return
    
    # Functon to run sheet_builder in 'sheet' mode
    def run_in_sheet_mode(self):
        self.run_sheet_builder('sheet')
        return
    
    # Functon to run sheet_builder in 'symmetric-tube' mode
    def run_in_symmetric_tube_mode(self):
        self.run_sheet_builder('symmetric-tube')
        return
    
    # Functon to run sheet_builder in 'chiral-tube' mode
    def run_in_chiral_tube_mode(self):
        self.run_sheet_builder('chiral-tube')
        return
    
    # Function to run imported sheet_builder main() functions
    def run_sheet_builder(self, run_mode): 
        # Set up log
        log = io_functions.LUNAR_logger()
        log.configure(level='production', print2console=True,  write2log=True)
        valid_inputs = True
        
        # Get information from GUI
        boolean = {'False':False, 'True':True}
        sheet_basename = self.sheet_basename.get()
        symmetric_tube_basename = self.symmetric_tube_basename.get()
        chiral_tube_basename = self.chiral_tube_basename.get()
        parent_directory = self.parent_directory.get() 
        try: length_in_perpendicular = float(self.length_in_perpendicular.get())
        except:
            length_in_perpendicular = 5
            log.GUI_error(f'ERROR perpendicular length is not a float {self.length_in_perpendicular.get()}')
            valid_inputs = False
        try: length_in_edgetype = float(self.length_in_edgetype.get())
        except:
            length_in_edgetype = 5
            log.GUI_error(f'ERROR perpendicular length is not a float {self.length_in_edgetype.get()}')
            valid_inputs = False
        sheet_edgetype = self.sheet_edgetype.get()
        types = {1:self.type1.get(),
                 2:self.type2.get(),
                 3:self.type3.get(),
                 4:self.type4.get()}
        try: bond_length = float(self.bond_length.get())
        except:
            bond_length = 1.4
            log.GUI_error(f'ERROR bond length is not a float {self.bond_length.get()}')
            valid_inputs = False
        try: sheet_layer_spacing = float(self.sheet_layer_spacing.get())
        except:
            sheet_layer_spacing = 3.4
            log.GUI_error(f'ERROR sheet layer spacing is not a float {self.sheet_layer_spacing.get()}')
            valid_inputs = False
        try: sheet_nlayers = int(self.sheet_nlayers.get())
        except:
            sheet_nlayers = 1
            log.GUI_error(f'ERROR sheet nlayers is not an int {self.sheet_nlayers.get()}')
            valid_inputs = False
        stacking = self.stacking.get()
        plane = self.plane.get()
        tube_edgetype = self.tube_edgetype.get()
        try: tube_layer_spacing = float(self.tube_layer_spacing.get())
        except:
            tube_layer_spacing = 3.4
            log.GUI_error(f'ERROR tube layer spacingis not a float {self.tube_layer_spacing.get()}')
            valid_inputs = False 
        try: symmetric_ntubes = int(self.symmetric_ntubes.get())
        except:
            symmetric_ntubes = 1
            log.GUI_error(f'ERROR symmetric_ is not an int {self.symmetric_ntubes.get()}')
            valid_inputs = False
        try: symmetric_length = float(self.symmetric_length.get())
        except:
            symmetric_length = 20
            log.GUI_error(f'ERROR sheet layer spacing is not a float {self.symmetric_length.get()}')
            valid_inputs = False 
        try: diameter = float(self.diameter.get())
        except:
            diameter = 12
            log.GUI_error(f'ERROR diameter is not a float {self.diameter.get()}')
            valid_inputs = False
        try: n = int(self.n.get())
        except:
            n = 5
            log.GUI_error(f'ERROR sheet n is not an int {self.n.get()}')
            valid_inputs = False  
        try: m = int(self.m.get())
        except:
            m = 5
            log.GUI_error(f'ERROR sheet m is not an int {self.m.get()}')
            valid_inputs = False 
        try: 
            chiral_length = float(self.chiral_length.get())
        except:
            chiral_length = 20
            log.GUI_error(f'ERROR chiral length is not a float {self.chiral_length.get()}')
            valid_inputs = False 
        symmetric_tube_axis = self.symmetric_tube_axis.get()
        chiral_tube_axis = self.chiral_tube_axis.get()
        find_bonds = boolean[self.find_bonds.get()]
        periodic_bonds = boolean[self.periodic_bonds.get()]
        charges = self.charges
        masses = self.masses
        cutter = [i.get() for i in self.cutter if i.get()]
        commandline_inputs = []
        
        try: seed = int(self.seed.get())
        except:
            seed = 0
            log.GUI_error(f'ERROR functional_seed is not an int {self.functional_seed.get()}')
            valid_inputs = False 
        functional_atoms = self.functional_atoms.get()
        terminating_atoms = self.terminating_atoms.get()
        grafting_files = self.grafting_files.get()
        try: minimum_distance = float(self.minimum_distance.get())
        except: minimum_distance = self.minimum_distance.get()
        
        # Run LUNAR/sheet_builder
        if valid_inputs:
            try:                
                inputs = (sheet_basename, symmetric_tube_basename, chiral_tube_basename, run_mode, parent_directory, length_in_perpendicular, length_in_edgetype, sheet_edgetype, types,
                          bond_length, sheet_layer_spacing, sheet_nlayers, stacking, plane, tube_edgetype, tube_layer_spacing, symmetric_ntubes, symmetric_length, diameter, n, m,
                          chiral_length, symmetric_tube_axis, chiral_tube_axis, find_bonds, periodic_bonds, charges, masses, seed, functional_atoms, terminating_atoms, 
                          grafting_files, minimum_distance, cutter, commandline_inputs, log)
                
                t1=threading.Thread(target=main, args=inputs)
                t1.start()
                t1.join()
            except Exception:
                log.GUI_error(traceback.format_exc())
        self.popup(log.logged, title='Outputs')
        return     
    
    # Function to pop-up scrollable text
    def popup(self, out, title='Outputs'):
        page = Toplevel(self.root)
        page.title(title)
        outputs = ScrolledText(page, height=30, width=150, font=('consolas', '12', 'normal'))
        outputs.pack()
        outputs.insert(tk.INSERT, '\n'.join(out))
        outputs.config(state=tk.DISABLED)
        return
    
    # Function to update py script default settings
    def update_py_script(self):
        # Get information from GUI
        boolean = {'False':False, 'True':True}
        sheet_basename = self.sheet_basename.get()
        symmetric_tube_basename = self.symmetric_tube_basename.get()
        chiral_tube_basename = self.chiral_tube_basename.get()
        parent_directory = io_functions.path_to_string(self.parent_directory.get())
        length_in_perpendicular = float(self.length_in_perpendicular.get())
        length_in_edgetype = float(self.length_in_edgetype.get())
        sheet_edgetype = self.sheet_edgetype.get()
        types = {1:self.type1.get(),
                 2:self.type2.get(),
                 3:self.type3.get(),
                 4:self.type4.get()}
        bond_length = float(self.bond_length.get())
        sheet_layer_spacing = float(self.sheet_layer_spacing.get())
        sheet_nlayers = int(self.sheet_nlayers.get())
        stacking = self.stacking.get()
        plane = self.plane.get()
        tube_edgetype = self.tube_edgetype.get()
        tube_layer_spacing = float(self.tube_layer_spacing.get())
        symmetric_ntubes = int(self.symmetric_ntubes.get())
        symmetric_length = float(self.symmetric_length.get())
        diameter = float(self.diameter.get())
        n = int(self.n.get()) 
        m = int(self.m.get())
        chiral_length = float(self.chiral_length.get())
        symmetric_tube_axis = self.symmetric_tube_axis.get()
        chiral_tube_axis = self.chiral_tube_axis.get()
        find_bonds = boolean[self.find_bonds.get()]
        periodic_bonds = boolean[self.periodic_bonds.get()]
        seed = int(self.seed.get())
        functional_atoms = self.functional_atoms.get()
        terminating_atoms = self.terminating_atoms.get()
        grafting_files = self.grafting_files.get()
        try: minimum_distance = float(self.minimum_distance.get())
        except: minimum_distance = self.minimum_distance.get()
        cutter = [i.get() for i in self.cutter if i.get()]
        
        # Read current py script and re-write with new settings
        print('Updating settings in: {}, from current GUI settings'.format(self.filename))
        lines = psm.read(self.filename)
        with open(self.filename, 'w') as f:
            inputsflag = True
            for line in lines:
                if line.startswith('sheet_basename') and inputsflag:
                    line = psm.parse_and_modify(line, sheet_basename, stringflag=True, splitchar='=')
                if line.startswith('symmetric_tube_basename') and inputsflag:
                    line = psm.parse_and_modify(line, symmetric_tube_basename, stringflag=True, splitchar='=')
                if line.startswith('chiral_tube_basename') and inputsflag:
                    line = psm.parse_and_modify(line, chiral_tube_basename, stringflag=True, splitchar='=')
                if line.startswith('parent_directory') and inputsflag:
                    line = psm.parse_and_modify(line, parent_directory, stringflag=True, splitchar='=')
                if line.startswith('length_in_perpendicular') and inputsflag:
                    line = psm.parse_and_modify(line, length_in_perpendicular, stringflag=False, splitchar='=')
                if line.startswith('length_in_edgetype ') and inputsflag:
                    line = psm.parse_and_modify(line, length_in_edgetype, stringflag=False, splitchar='=')
                if line.startswith('sheet_edgetype') and inputsflag:
                    line = psm.parse_and_modify(line, sheet_edgetype, stringflag=True, splitchar='=')
                if line.startswith('types') and inputsflag:
                    line = psm.parse_and_modify(line, str(types), stringflag=False, splitchar='=')
                    
                    
                if line.startswith('cutter') and inputsflag:
                    line = psm.parse_and_modify(line, str(cutter), stringflag=False, splitchar='=')
                    
                if line.startswith('bond_length') and inputsflag:
                    line = psm.parse_and_modify(line, bond_length, stringflag=False, splitchar='=')
                if line.startswith('sheet_layer_spacing') and inputsflag:
                    line = psm.parse_and_modify(line, sheet_layer_spacing, stringflag=False, splitchar='=')
                if line.startswith('sheet_nlayers') and inputsflag:
                    line = psm.parse_and_modify(line, sheet_nlayers, stringflag=False, splitchar='=')
                if line.startswith('stacking') and inputsflag:
                    line = psm.parse_and_modify(line, stacking, stringflag=True, splitchar='=')
                if line.startswith('plane') and inputsflag:
                    line = psm.parse_and_modify(line, plane, stringflag=True, splitchar='=')
                if line.startswith('tube_edgetype') and inputsflag:
                    line = psm.parse_and_modify(line, tube_edgetype, stringflag=True, splitchar='=')
                if line.startswith('tube_layer_spacing') and inputsflag:
                    line = psm.parse_and_modify(line, tube_layer_spacing, stringflag=False, splitchar='=') 
                if line.startswith('symmetric_ntubes') and inputsflag:
                    line = psm.parse_and_modify(line, symmetric_ntubes, stringflag=False, splitchar='=')
                if line.startswith('symmetric_length') and inputsflag:
                    line = psm.parse_and_modify(line, symmetric_length, stringflag=False, splitchar='=')  
                if line.startswith('diameter') and inputsflag:
                    line = psm.parse_and_modify(line, diameter, stringflag=False, splitchar='=')   
                if line.startswith('n') and inputsflag:
                    line = psm.parse_and_modify(line, n, stringflag=False, splitchar='=')   
                if line.startswith('m') and inputsflag and 'masses' not in line:
                    line = psm.parse_and_modify(line, m, stringflag=False, splitchar='=')
                if line.startswith('chiral_length') and inputsflag:
                    line = psm.parse_and_modify(line, chiral_length, stringflag=False, splitchar='=')
                if line.startswith('symmetric_tube_axis') and inputsflag:
                    line = psm.parse_and_modify(line, symmetric_tube_axis, stringflag=True, splitchar='=')    
                if line.startswith('chiral_tube_axis') and inputsflag:
                    line = psm.parse_and_modify(line, chiral_tube_axis, stringflag=True, splitchar='=')
                if line.startswith('find_bonds') and inputsflag:
                    line = psm.parse_and_modify(line, find_bonds, stringflag=False, splitchar='=')
                if line.startswith('periodic_bonds') and inputsflag:
                    line = psm.parse_and_modify(line, periodic_bonds, stringflag=False, splitchar='=')
                if line.startswith('seed') and inputsflag:
                    line = psm.parse_and_modify(line, seed, stringflag=False, splitchar='=')
                if line.startswith('functional_atoms') and inputsflag:
                    line = psm.parse_and_modify(line, functional_atoms, stringflag=True, splitchar='=')
                if line.startswith('terminating_atoms') and inputsflag:
                    line = psm.parse_and_modify(line, terminating_atoms, stringflag=True, splitchar='=')
                if line.startswith('grafting_files') and inputsflag:
                    line = psm.parse_and_modify(line, grafting_files, stringflag=True, splitchar='=')
                if line.startswith('minimum_distance') and inputsflag and isinstance(minimum_distance, (float, int)):
                    line = psm.parse_and_modify(line, minimum_distance, stringflag=False, splitchar='=') 
                if line.startswith('minimum_distance') and inputsflag and isinstance(minimum_distance, str):
                    line = psm.parse_and_modify(line, minimum_distance, stringflag=True, splitchar='=') 
                if line.startswith('if __name__ == "__main__":'): inputsflag = False
                f.write(line)
        return