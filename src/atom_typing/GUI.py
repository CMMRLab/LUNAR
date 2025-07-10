# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
May 13th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.GUI_scale_settings as GUI_scale_settings
import src.io_functions as io_functions
from src.atom_typing.main import main
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


#########################
# LUNAR/atom_typing GUI #
#########################
class atom_typing_GUI:
    def __init__(self, topofile, bondfile, parent_directory, newfile, ff_name, delete_atoms, mass_map, bondorder, maxbonded, boundary,
                       vdw_radius_scale, reset_charges, print_options, commandline_inputs, bonds_via_distance_override, pdb_file,
                       chargefile, include_comments_nta, GUI_zoom):
        
        # Pass certain inputs as attribute (These will only be able to be adjusted from the python script)
        self.mass_map = mass_map
        self.bondorder = bondorder
        self.maxbonded = maxbonded
        
        # Find present working directory
        self.pwd = os.getcwd()
        self.filepath = self.pwd
        self.filename = os.path.join(self.pwd, 'atom_typing.py')

        # Initialize window
        self.root = tk.Tk()
        self.root.title('LUNAR/atom_typing.py GUI v1.0')
        self.root.resizable(width=False, height=False)
        #self.root.geometry('600x400')
        
        # Initalize main frame
        self.frame = tk.Frame(self.root)
        self.frame.pack()
        
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
        self.maxwidth = 100

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
        font_settings = (self.font_type, font_size)
        
        #--------------#
        # Inputs frame #
        #--------------#
        # Initalize  inputs frame
        self.inputs_frame = tk.LabelFrame(self.frame, text='Inputs', font=font_settings)
        self.inputs_frame.grid(row=0, column=0, columnspan=2, padx=xpadding, pady=ypadding)
        
        # topofile selection button
        self.topofile = tk.Entry(self.inputs_frame, width=int(1.2*maxwidth), font=font_settings)
        self.topofile.insert(0, topofile)
        self.topofile.grid(column=1, row=0)
        self.topofile_button = tk.Button(self.inputs_frame, text='topofile', font=font_settings, command=self.topofile_path)
        self.topofile_button.grid(column=0, row=0)
        
        # bondfile selection button
        self.bondfile = tk.Entry(self.inputs_frame, width=int(1.2*maxwidth), font=font_settings)
        self.bondfile.insert(0, bondfile)
        self.bondfile.grid(column=1, row=1)
        self.bondfile_button = tk.Button(self.inputs_frame, text='bondfile', font=font_settings, command=self.bondfile_path)
        self.bondfile_button.grid(column=0, row=1)
        
        # chargefile selection button
        self.chargefile = tk.Entry(self.inputs_frame, width=int(1.2*maxwidth), font=font_settings)
        self.chargefile.insert(0, chargefile)
        self.chargefile.grid(column=1, row=2)
        self.chargefile_button = tk.Button(self.inputs_frame, text='chargefile', font=font_settings, command=self.chargefile_path)
        self.chargefile_button.grid(column=0, row=2)
        
        # parent_directory entry
        self.parent_directory = tk.Entry(self.inputs_frame, width=int(1.2*maxwidth), font=font_settings)
        self.parent_directory.insert(0, parent_directory)
        self.parent_directory.grid(column=1, row=3)
        self.dir_button = tk.Button(self.inputs_frame, text='parent_directory', font=font_settings, command=self.directory_path)
        self.dir_button.grid(column=0, row=3)
        
        # newfile entry
        self.newfile = tk.Entry(self.inputs_frame, width=int(1.2*maxwidth), font=font_settings)
        self.newfile.insert(0, newfile)
        self.newfile.grid(column=1, row=5)
        self.newfile_label = tk.Label(self.inputs_frame, text='newfile', font=font_settings)
        self.newfile_label.grid(column=0, row=5)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.inputs_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=ypadding)
            
        
        #---------------#
        # Options frame #
        #---------------#
        # Initalize  options frame
        self.options_frame = tk.LabelFrame(self.frame, text='Options', font=font_settings)
        self.options_frame.grid(row=1, column=0, columnspan=2, sticky='news', padx=xpadding, pady=ypadding)
                
        # ff_class drop down menu
        styles = ['PCFF-IFF', 'PCFF', 'compass', 'CVFF-IFF', 'CVFF', 'Clay-FF', 'DREIDING', 'OPLS-AA', 'general:0',
                  'general:1', 'general:2', 'general:3', 'general:4', 'general-pp:2', 'general-pp:3', 'general-pp:4']
        self.ff_name = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/7), font=font_settings)
        self.ff_name.current(styles.index(ff_name))
        self.ff_name.grid(column=0, row=1)
        self.ff_name_label = tk.Label(self.options_frame, text='ff_name', font=font_settings)
        self.ff_name_label.grid(column=0, row=0)
        
        # reset_charges drop down menu
        styles = [True, False]
        self.reset_charges = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/10), font=font_settings)
        self.reset_charges.current(styles.index(reset_charges))
        self.reset_charges.grid(column=1, row=1)
        self.reset_charges_label = tk.Label(self.options_frame, text='reset_charges', font=font_settings)
        self.reset_charges_label.grid(column=1, row=0)
        
        # delete_atoms method drop down menu
        method = ['mass', 'size']
        self.da_method = ttk.Combobox(self.options_frame, values=method, width=int(maxwidth/7.5), font=font_settings)
        self.da_method.current(method.index(delete_atoms['method']))
        self.da_method.grid(column=2, row=1)
        self.da_method_label = tk.Label(self.options_frame, text="delete_atoms['method']", font=font_settings)
        self.da_method_label.grid(column=2, row=0)
        
        # delete_atoms criteria entry
        self.da_criteria = tk.Entry(self.options_frame, width=int(maxwidth/7.5), font=font_settings)
        self.da_criteria.insert(0, delete_atoms['criteria'])
        self.da_criteria.grid(column=3, row=1)
        self.da_criteria_label = tk.Label(self.options_frame, text="delete_atoms['criteria']", font=font_settings)
        self.da_criteria_label.grid(column=3, row=0)
        
        # pdb_file method drop down menu
        method = ['skip', 'types', 'typeIDs']
        self.pdb_file = ttk.Combobox(self.options_frame, values=method, width=int(maxwidth/12), font=font_settings)
        self.pdb_file.current(method.index(pdb_file))
        self.pdb_file.grid(column=4, row=1)
        self.pdb_file_label = tk.Label(self.options_frame, text='pdb_file', font=font_settings)
        self.pdb_file_label.grid(column=4, row=0)

        # include_comments_nta Boolean drop down menu
        styles = [True, False]
        self.include_comments_nta = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/7.5), font=font_settings)
        self.include_comments_nta.current(styles.index(include_comments_nta))
        self.include_comments_nta.grid(column=5, row=1)
        self.include_comments_nta_label = tk.Label(self.options_frame, text='include_comments_nta', font=font_settings)
        self.include_comments_nta_label.grid(column=5, row=0)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.options_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))
            
            
        #--------------------#
        # Bond finding frame #
        #--------------------#
        # Initalize  options frame
        self.bond_frame = tk.LabelFrame(self.frame, text='Bond finding options', font=font_settings)
        self.bond_frame.grid(row=2, column=0, columnspan=2, sticky='news', padx=xpadding, pady=ypadding)
        
        # vdw_radius_scale entry
        self.vdw_radius_scale = tk.Entry(self.bond_frame, width=int(maxwidth/7.5), font=font_settings)
        self.vdw_radius_scale.insert(0, vdw_radius_scale)
        self.vdw_radius_scale.grid(column=1, row=1)
        self.vdw_radius_scale_label = tk.Label(self.bond_frame, text='vdw_radius_scale', font=font_settings)
        self.vdw_radius_scale_label.grid(column=1, row=0)
        
        # boundary entry
        self.boundary = tk.Entry(self.bond_frame, width=int(maxwidth/7.5), font=font_settings)
        self.boundary.insert(0, boundary)
        self.boundary.grid(column=2, row=1)
        self.boundary_label = tk.Label(self.bond_frame, text='boundary', font=font_settings)
        self.boundary_label.grid(column=2, row=0)
        
        # bonds_via_distance_override drop down menu
        styles = [True, False]
        self.bonds_via_distance_override = ttk.Combobox(self.bond_frame, values=styles, width=int(maxwidth/7.5), font=font_settings)
        self.bonds_via_distance_override.current(styles.index(bonds_via_distance_override))
        self.bonds_via_distance_override.grid(column=3, row=1)
        self.bonds_via_distance_override_label = tk.Label(self.bond_frame, text='bonds_via_distance_override', font=font_settings)
        self.bonds_via_distance_override_label.grid(column=3, row=0)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.bond_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=5)
        
        #------------#
        # Run button #
        #------------#
        self.run = tk.Button(self.frame, text='Run LUNAR/atom_typing.py', font=font_settings, command=self.run_LUNAR)
        self.run.grid(row=3, column=0, columnspan=2, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
        #-----------------#
        # update defaults #
        #-----------------#
        self.update = tk.Button(self.frame, text='Save the current GUI settings as the default GUI settings', font=font_settings, command=self.update_py_script)
        self.update.grid(row=4, column=0, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
        #------------#
        # Quick help #
        #------------#
        self.quick_help = tk.Button(self.frame, text='Quick help', font=font_settings, command=self.quickhelp)
        self.quick_help.grid(row=4, column=1, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
        
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
            txt = os.path.join(self.pwd, 'src/GUI_quick_help_pages/atom_typing.txt')
            logged = []
            with open(txt, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    # Strip comment's and split by whitespace
                    line = line.split('#')[0]
                    line = line.rstrip()
                    logged.append(line)
        except: # except something failed
            logged.append('FAILED to read LUNAR/src/GUI_quick_help_pages/atom_typing.txt document.')
            logged.append('Most likely cause is the atom_typing.txt file was renamed in the')
            logged.append('LUNAR/src/GUI_quick_help_pages/ directory or directory names were changed.')
        self.popup(logged, title='Quick help')
        return
    
    # Function to get filepath for topofile
    def topofile_path(self):
        path = ''
        try:
            ftypes = (('all files', '*.*'), ('data files', '*.data *.data.gz'), ('mol files', '*.mol'),
                      ('mol2 files', '*.mol2'), ('mdf files', '*.mdf'), ('sdf files', '*.sdf'),
                      ('pdb files', '*.pdb'))
            path = filedialog.askopenfilename(title='Open topofile?', filetypes=ftypes)
        except: path = filedialog.askopenfilename(title='Open topofile?')
        if path:
            self.filepath = os.path.dirname(os.path.abspath(path))
            path = os.path.relpath(path)
            self.topofile.delete(0, tk.END); self.topofile.insert(0, path);
        return
    
    # Function to get filepath for bondfile
    def bondfile_path(self):
        path = filedialog.askopenfilename(title='Open bondfile?')
        if path:
            self.filepath = os.path.dirname(os.path.abspath(path))
            path = os.path.relpath(path)
            self.bondfile.delete(0, tk.END); self.bondfile.insert(0, path);
        return
    
    # Function to get filepath for chargefile
    def chargefile_path(self):
        ftypes = (('txt files', '*.txt'), ('all files', '*.*'))
        path = filedialog.askopenfilename(title='Open chargefile?', filetypes=ftypes)
        if path:
            self.filepath = os.path.dirname(os.path.abspath(path))
            path = os.path.relpath(path)
            self.chargefile.delete(0, tk.END); self.chargefile.insert(0, path);
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
        print('Terminating atom_typing GUI'); self.root.destroy();
        return
    
    # Function to run imported LUNAR functions
    def run_LUNAR(self):  
        # Set up log
        log = io_functions.LUNAR_logger()
        valid_inputs = True
        
        # Get information from GUI
        boolean = {'False':False, 'True':True}
        topofile = self.topofile.get()
        bondfile = self.bondfile.get()
        chargefile = self.chargefile.get()
        parent_directory = self.parent_directory.get() 
        newfile = self.newfile.get()
        ff_name = self.ff_name.get()
        reset_charges = boolean[self.reset_charges.get()]       
        delete_atoms = {'method': self.da_method.get(),
                        'criteria': float(self.da_criteria.get()) }
        vdw_radius_scale = float(self.vdw_radius_scale.get())
        boundary = self.boundary.get()
        bonds_via_distance_override = boolean[self.bonds_via_distance_override.get()]
        print_options = False
        mass_map = self.mass_map
        bondorder = self.bondorder
        maxbonded = self.maxbonded
        commandline_inputs = []
        pdb_file = self.pdb_file.get()
        include_comments_nta = boolean[self.include_comments_nta.get()]

        # Run LUNAR/atom_typing
        if valid_inputs:
            try: 
                inputs = (topofile, bondfile, parent_directory, newfile, ff_name, delete_atoms, mass_map, bondorder, maxbonded, boundary,
                          vdw_radius_scale, reset_charges, print_options, commandline_inputs, bonds_via_distance_override, pdb_file, chargefile,
                          include_comments_nta, log)            
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
        page.resizable(width=False, height=False)
        return


    # Function to update py script default settings
    def update_py_script(self):
        # Get information from GUI
        boolean = {'False':False, 'True':True}
        topofile = io_functions.path_to_string(self.topofile.get())
        bondfile = io_functions.path_to_string(self.bondfile.get())
        chargefile = io_functions.path_to_string(self.chargefile.get())
        parent_directory = io_functions.path_to_string(self.parent_directory.get()) 
        newfile = self.newfile.get()
        ff_name = self.ff_name.get()
        reset_charges = boolean[self.reset_charges.get()]       
        vdw_radius_scale = float(self.vdw_radius_scale.get())
        boundary = self.boundary.get()
        bonds_via_distance_override = boolean[self.bonds_via_distance_override.get()]
        delete_atoms = {'method': self.da_method.get(),
                        'criteria': float(self.da_criteria.get()) }
        pdb_file = self.pdb_file.get()
        include_comments_nta = boolean[self.include_comments_nta.get()]
        
        # Read current py script and re-write with new settings
        print('Updating settings in: {}, from current GUI settings'.format(self.filename))
        lines = psm.read(self.filename)
        with open(self.filename, 'w') as f:
            inputsflag = True; da_flag = False
            for line in lines:
                # if line.startswith('use_GUI') and inputsflag:
                #     line = psm.parse_and_modify(line, True, stringflag=False, splitchar='=')
                # if line.startswith('topofile') and inputsflag:
                #     line = psm.parse_and_modify(line, topofile, stringflag=True, splitchar='=')
                # if line.startswith('bondfile') and inputsflag:
                #     line = psm.parse_and_modify(line, bondfile, stringflag=True, splitchar='=')
                if line.startswith('chargefile') and inputsflag:
                    line = psm.parse_and_modify(line, chargefile, stringflag=True, splitchar='=')
                if line.startswith('parent_directory') and inputsflag:
                    line = psm.parse_and_modify(line, parent_directory, stringflag=True, splitchar='=')
                if line.startswith('pdb_file') and inputsflag:
                    line = psm.parse_and_modify(line, pdb_file, stringflag=True, splitchar='=')            
                if line.startswith('newfile') and inputsflag:
                    line = psm.parse_and_modify(line, newfile, stringflag=True, splitchar='=')    
                if line.startswith('ff_name') and inputsflag:
                    line = psm.parse_and_modify(line, ff_name, stringflag=True, splitchar='=')
                if line.startswith('reset_charges') and inputsflag:
                    line = psm.parse_and_modify(line, reset_charges, stringflag=False, splitchar='=')
                if line.startswith('include_comments_nta') and inputsflag:
                    line = psm.parse_and_modify(line, include_comments_nta, stringflag=False, splitchar='=')
                if line.startswith('vdw_radius_scale') and inputsflag:
                    line = psm.parse_and_modify(line, vdw_radius_scale, stringflag=False, splitchar='=')
                if line.startswith('boundary') and inputsflag:
                    line = psm.parse_and_modify(line, boundary, stringflag=True, splitchar='=')  
                if line.startswith('bonds_via_distance_override') and inputsflag:
                    line = psm.parse_and_modify(line, bonds_via_distance_override, stringflag=False, splitchar='=')
                if line.startswith('delete_atoms') and inputsflag: da_flag = True
                if 'method' in line and inputsflag and da_flag and not line.startswith('#'):
                    line = psm.parse_and_modify(line, delete_atoms['method'], stringflag=True, splitchar=':')
                if 'criteria' in line and inputsflag and da_flag and not line.startswith('#'):
                    line = psm.parse_and_modify(line, delete_atoms['criteria'], stringflag=False, splitchar=':')
                if line.startswith('if __name__ == "__main__":'): inputsflag = False
                f.write(line)
        return    