# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
November 13th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.io_functions as io_functions
import src.py_script_modifier as psm
from src.all2lmp.main import main
from tkinter.scrolledtext import ScrolledText
from tkinter import filedialog
from tkinter import Toplevel
from tkinter import ttk
import tkinter as tk
import threading
import traceback
import math
import os




#####################
# LUNAR/all2lmp GUI #
#####################
class all2lmp_GUI:
    def __init__(self, topofile, nta_file, frc_file, assumed, parent_directory, newfile, atom_style, ff_class, use_auto_equivalence, use_assumed_auto_fill, reset_molids, 
                 reset_charges, write_txt_comments, write_bond_react, print_options, use_morse_bonds, include_type_labels, add2box, ignore_missing_parameters, 
                 shift, rotate, GUI_zoom, commandline_inputs):
        
        # Find present working directory
        self.pwd = os.getcwd()
        self.filepath = self.pwd
        self.filename = os.path.join(self.pwd, 'all2lmp.py')

        # Initialize window
        self.root = tk.Tk()
        self.root.title('LUNAR/all2lmp.py GUI v1.0')
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
        self.xpadding = 20
        self.ypadding = 10
        self.maxwidth = 120
        
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
        self.topofile = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.topofile.insert(0, topofile)
        self.topofile.grid(column=1, row=0)
        self.topofile_button = tk.Button(self.inputs_frame, text='topofile', font=font_settings, command=self.topofile_path)
        self.topofile_button.grid(column=0, row=0)
        
        # nta_file selection button
        self.nta_file = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.nta_file.insert(0, nta_file)
        self.nta_file.grid(column=1, row=1)
        self.nta_file_button = tk.Button(self.inputs_frame, text='nta_file', font=font_settings, command=self.nta_file_path)
        self.nta_file_button.grid(column=0, row=1)
        
        # frc_file selection button
        self.frc_file = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.frc_file.insert(0, frc_file)
        self.frc_file.grid(column=1, row=2)
        self.frc_file_button = tk.Button(self.inputs_frame, text='frc_file', font=font_settings, command=self.frc_file_path)
        self.frc_file_button.grid(column=0, row=2)
        
        # assumed selection button
        self.assumed = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.assumed.insert(0, assumed)
        self.assumed.grid(column=1, row=3)
        self.assumed_button = tk.Button(self.inputs_frame, text='assumed', font=font_settings, command=self.assumed_path)
        self.assumed_button.grid(column=0, row=3)
        
        # parent_directory entry
        self.parent_directory = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.parent_directory.insert(0, parent_directory)
        self.parent_directory.grid(column=1, row=4)
        self.dir_button = tk.Button(self.inputs_frame, text='parent_directory', font=font_settings, command=self.directory_path)
        self.dir_button.grid(column=0, row=4)
        
        # newfile entry
        self.newfile = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
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
        styles = ['0', '1', '2', 'i', 'ilmp', 'd', 's1', 's2']
        self.ff_class = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/12), font=font_settings)
        self.ff_class.current(styles.index(ff_class))
        self.ff_class.grid(column=0, row=1)
        self.ff_class_label = tk.Label(self.options_frame, text='ff_class', font=font_settings)
        self.ff_class_label.grid(column=0, row=0)
        
        # atom_style drop down
        styles = ['full', 'charge', 'molecular', 'angle', 'bond', 'atomic', 'dipole', 'dpd', 'line']
        self.atom_style = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/12), font=font_settings)
        self.atom_style.current(styles.index(atom_style))
        self.atom_style.grid(column=1, row=1)
        self.atom_style_label = tk.Label(self.options_frame, text='atom_style', font=font_settings)
        self.atom_style_label.grid(column=1, row=0)
        
        # use_auto_equivalence drop down menu
        styles = [True, False]
        self.use_auto_equivalence = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/12), font=font_settings)
        self.use_auto_equivalence.current(styles.index(use_auto_equivalence))
        self.use_auto_equivalence.grid(column=2, row=1)
        self.use_auto_equivalence_label = tk.Label(self.options_frame, text='use_auto_equivalence', font=font_settings)
        self.use_auto_equivalence_label.grid(column=2, row=0)
        
        # use_morse_bonds drop down menu
        styles = [True, False]
        self.use_morse_bonds = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/12), font=font_settings)
        self.use_morse_bonds.current(styles.index(use_morse_bonds))
        self.use_morse_bonds.grid(column=3, row=1)
        self.use_morse_bonds_label = tk.Label(self.options_frame, text='use_morse_bonds', font=font_settings)
        self.use_morse_bonds_label.grid(column=3, row=0)
        
        # use_assumed_auto_fill drop down menu
        styles = [True, False]
        self.use_assumed_auto_fill = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/12), font=font_settings)
        self.use_assumed_auto_fill.current(styles.index(use_assumed_auto_fill))
        self.use_assumed_auto_fill.grid(column=4, row=1)
        self.use_assumed_auto_fill_label = tk.Label(self.options_frame, text='use_assumed_auto_fill', font=font_settings)
        self.use_assumed_auto_fill_label.grid(column=4, row=0)
        
        # fav entry
        self.add2box = tk.Entry(self.options_frame, width=int(maxwidth/12), font=font_settings)
        self.add2box.insert(0, add2box)
        self.add2box.grid(column=5, row=1)
        self.add2box_label = tk.Label(self.options_frame, text='add2box (float or int)', font=font_settings)
        self.add2box_label.grid(column=5, row=0)
        
        # reset_molids drop down menu
        styles = [True, False]
        self.reset_molids = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/12), font=font_settings)
        self.reset_molids.current(styles.index(reset_molids))
        self.reset_molids.grid(column=0, row=3)
        self.reset_molids_label = tk.Label(self.options_frame, text='reset_molids', font=font_settings)
        self.reset_molids_label.grid(column=0, row=2)
        
        # reset_charges drop down menu
        styles = [True, False]
        self.reset_charges = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/12), font=font_settings)
        self.reset_charges.current(styles.index(reset_charges))
        self.reset_charges.grid(column=1, row=3)
        self.reset_charges_label = tk.Label(self.options_frame, text='reset_charges', font=font_settings)
        self.reset_charges_label.grid(column=1, row=2)
        
        # write_txt_comments drop down menu
        styles = [True, False]
        self.write_txt_comments = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/12), font=font_settings)
        self.write_txt_comments.current(styles.index(write_txt_comments))
        self.write_txt_comments.grid(column=2, row=3)
        self.write_txt_comments_label = tk.Label(self.options_frame, text='write_txt_comments', font=font_settings)
        self.write_txt_comments_label.grid(column=2, row=2)
        
        # write_bond_react drop down menu
        styles = [True, False]
        self.write_bond_react = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/12), font=font_settings)
        self.write_bond_react.current(styles.index(write_bond_react))
        self.write_bond_react.grid(column=3, row=3)
        self.write_bond_react_label = tk.Label(self.options_frame, text='write_bond_react', font=font_settings)
        self.write_bond_react_label.grid(column=3, row=2)
        
        # include_type_labels drop down menu
        styles = [True, False]
        self.include_type_labels = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/12), font=font_settings)
        self.include_type_labels.current(styles.index(include_type_labels))
        self.include_type_labels.grid(column=4, row=3)
        self.include_type_labels_label = tk.Label(self.options_frame, text='include_type_labels', font=font_settings)
        self.include_type_labels_label.grid(column=4, row=2)
        
        # include_type_labels drop down menu
        styles = [True, False]
        self.ignore_missing_parameters = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/12), font=font_settings)
        self.ignore_missing_parameters.current(styles.index(ignore_missing_parameters))
        self.ignore_missing_parameters.grid(column=5, row=3)
        self.ignore_missing_parameters_label = tk.Label(self.options_frame, text='ignore_missing_parameters', font=font_settings)
        self.ignore_missing_parameters_label.grid(column=5, row=2)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.options_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))
        
        #----------------------#
        # Atom positions frame #
        #----------------------#
        # Initalize  options frame
        self.atom_positions_frame = tk.LabelFrame(self.frame, text='Atom positions', font=font_settings)
        self.atom_positions_frame.grid(row=2, column=0, columnspan=2, sticky='news', padx=xpadding, pady=ypadding)
        
        # X-shift
        self.xs = ttk.Entry(self.atom_positions_frame, width=int(maxwidth/7), font=font_settings)
        self.xs.insert(0, shift['x'])
        self.xs.grid(column=0, row=1)
        self.xs_label = tk.Label(self.atom_positions_frame, text='X-shift', font=font_settings)
        self.xs_label.grid(column=0, row=0)
        
        # Y-shift
        self.ys = ttk.Entry(self.atom_positions_frame, width=int(maxwidth/7), font=font_settings)
        self.ys.insert(0, shift['y'])
        self.ys.grid(column=1, row=1)
        self.ys_label = tk.Label(self.atom_positions_frame, text='Y-shift', font=font_settings)
        self.ys_label.grid(column=1, row=0)
        
        # Z-shift
        self.zs = ttk.Entry(self.atom_positions_frame, width=int(maxwidth/7), font=font_settings)
        self.zs.insert(0, shift['z'])
        self.zs.grid(column=2, row=1)
        self.zs_label = tk.Label(self.atom_positions_frame, text='Z-shift', font=font_settings)
        self.zs_label.grid(column=2, row=0)
        
        # X-rotate
        self.xr = ttk.Entry(self.atom_positions_frame, width=int(maxwidth/7), font=font_settings)
        self.xr.insert(0, rotate['x'])
        self.xr.grid(column=3, row=1)
        self.xr_label = tk.Label(self.atom_positions_frame, text='rotate about X', font=font_settings)
        self.xr_label.grid(column=3, row=0)
        
        # Y-rotate
        self.yr = ttk.Entry(self.atom_positions_frame, width=int(maxwidth/7), font=font_settings)
        self.yr.insert(0, rotate['y'])
        self.yr.grid(column=4, row=1)
        self.yr_label = tk.Label(self.atom_positions_frame, text='rotate about Y', font=font_settings)
        self.yr_label.grid(column=4, row=0)
        
        # Z-rotate
        self.zr = ttk.Entry(self.atom_positions_frame, width=int(maxwidth/7), font=font_settings)
        self.zr.insert(0, rotate['z'])
        self.zr.grid(column=5, row=1)
        self.zr_label = tk.Label(self.atom_positions_frame, text='rotate about Z', font=font_settings)
        self.zr_label.grid(column=5, row=0)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.atom_positions_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))
        
        #------------#
        # Run button #
        #------------#
        self.run = tk.Button(self.frame, text='Run LUNAR/all2lmp.py', font=font_settings, command=self.run_LUNAR)
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
            txt = os.path.join(self.pwd, 'src/GUI_quick_help_pages/all2lmp.txt')
            logged = []
            with open(txt, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    # Strip comment's and split by whitespace
                    line = line.split('#')[0]
                    line = line.rstrip()
                    logged.append(line)
        except: # except something failed
            logged.append('FAILED to read LUNAR/src/GUI_quick_help_pages/all2lmp.txt document.')
            logged.append('Most likely cause is the all2lmp.txt file was renamed in the')
            logged.append('LUNAR/src/GUI_quick_help_pages/ directory or directory names were changed.')
        self.popup(logged, title='Quick help')
        return
    
    # Function to get filepath for topofile
    def topofile_path(self):
        ftypes = (('all files', '*.*'), ('data files', '*.data *.data.gz'), ('mol files', '*.mol'),
                  ('mol2 files', '*.mol2'), ('mdf files', '*.mdf'), ('sdf files', '*.sdf'),
                  ('pdb files', '*.pdb'))
        path = filedialog.askopenfilename(initialdir=self.filepath, title='Open topofile?', filetypes=ftypes)
        if path:
            self.filepath = os.path.dirname(os.path.abspath(path))
            path = os.path.relpath(path)
            self.topofile.delete(0, tk.END); self.topofile.insert(0, path);
        return
    
    # Function to get directory
    def directory_path(self):
        path = filedialog.askdirectory(initialdir=self.pwd)
        if path:
            path = os.path.relpath(path)
            self.parent_directory.delete(0, tk.END); self.parent_directory.insert(0, path);
        return
    
    # Function to get filepath for nta_file
    def nta_file_path(self):
        ftypes = (('all files', '*.*'), ('nta files', '*.nta'), ('car files', '*.car'))
        path = filedialog.askopenfilename(initialdir=self.filepath, title='Open nta_file?', filetypes=ftypes)
        self.filepath = os.path.dirname(os.path.abspath(path))
        if path:
            path = os.path.relpath(path)
            self.nta_file.delete(0, tk.END); self.nta_file.insert(0, path);
        return
    
    # Function to get filepath for frc_file
    def frc_file_path(self):
        ftypes = (('frc files', '*.frc'), ('all files', '*.*'))
        path = filedialog.askopenfilename(initialdir=self.pwd, title='Open frc file?', filetypes=ftypes)
        if path:
            path = os.path.relpath(path)
            self.frc_file.delete(0, tk.END); self.frc_file.insert(0, path);
        return
    
    # Function to get filepath for assumed file
    def assumed_path(self):
        ftypes = (('all files', '*.*'), ('assumed files', '*.coeff'))
        path = filedialog.askopenfilename(initialdir=self.pwd, title='Open assumed coeffs file?', filetypes=ftypes)
        if path:
            path = os.path.relpath(path)
            self.assumed.delete(0, tk.END); self.assumed.insert(0, path);
        return
    
    # Closing command    
    def closing(self):
        print('Terminating all2lmp GUI'); self.root.destroy();
        return
    
    # Function to run imported LUNAR functions
    def run_LUNAR(self):  
        # Set up log
        log = io_functions.LUNAR_logger()
        valid_inputs = True
        
        # Get information from GUI
        boolean = {'False':False, 'True':True}
        topofile = self.topofile.get()
        nta_file = self.nta_file.get()
        frc_file = self.frc_file.get()
        assumed = self.assumed.get()
        parent_directory = self.parent_directory.get() 
        atom_style = self.atom_style.get()
        newfile = self.newfile.get()
        ff_class = self.ff_class.get()
        use_auto_equivalence = boolean[self.use_auto_equivalence.get()]
        use_morse_bonds = boolean[self.use_morse_bonds.get()]
        use_assumed_auto_fill = boolean[self.use_assumed_auto_fill.get()]
        try: reset_molids = boolean[self.reset_molids.get()]
        except: reset_molids = self.reset_molids.get()
        reset_charges = boolean[self.reset_charges.get()]
        write_txt_comments = boolean[self.write_txt_comments.get()]
        write_bond_react = boolean[self.write_bond_react.get()]
        include_type_labels = boolean[self.include_type_labels.get()]
        ignore_missing_parameters = boolean[self.ignore_missing_parameters.get()]
        add2box = float(self.add2box.get())
        commandline_inputs = []
        print_options = False
        try:
            shift = {'x': float(self.xs.get()),
                     'y': float(self.ys.get()),
                     'z': float(self.zs.get())}
        except:
            log.GUI_error('ERROR shift X, Y, or Z rotation is not a float')
            valid_inputs = False
        try:
            try: rx = float(self.xr.get())
            except: rx = self.xr.get()
            try: ry = float(self.yr.get())
            except: ry = self.yr.get()
            try: rz = float(self.zr.get())
            except: rz = self.zr.get()
            rotate = {'x': rx,
                      'y': ry,
                      'z': rz}
        except:
            log.GUI_error('ERROR max X, Y, or Z rotation is not a float or string')
            valid_inputs = False


        # Run LUNAR/all2lmp
        if valid_inputs:
            try: 
                inputs = (topofile, nta_file, frc_file, assumed, parent_directory, newfile, atom_style, ff_class,
                          use_auto_equivalence, use_assumed_auto_fill, reset_molids, reset_charges, write_txt_comments,
                          write_bond_react, print_options, use_morse_bonds, include_type_labels, add2box, ignore_missing_parameters,
                          shift, rotate, commandline_inputs, log)
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
        nta_file = io_functions.path_to_string(self.nta_file.get())
        frc_file = io_functions.path_to_string(self.frc_file.get())
        assumed = io_functions.path_to_string(self.assumed.get())
        parent_directory = io_functions.path_to_string(self.parent_directory.get()) 
        atom_style = self.atom_style.get()
        newfile = self.newfile.get()
        ff_class = self.ff_class.get()
        
        use_auto_equivalence = boolean[self.use_auto_equivalence.get()]
        use_morse_bonds = boolean[self.use_morse_bonds.get()]
        use_assumed_auto_fill = boolean[self.use_assumed_auto_fill.get()]
        try: reset_molids = boolean[self.reset_molids.get()]
        except: reset_molids = self.reset_molids.get()
        reset_charges = boolean[self.reset_charges.get()]
        write_txt_comments = boolean[self.write_txt_comments.get()]
        write_bond_react = boolean[self.write_bond_react.get()]
        include_type_labels = boolean[self.include_type_labels.get()]
        ignore_missing_parameters = boolean[self.ignore_missing_parameters.get()]
        add2box = float(self.add2box.get())
        shift = {'x': float(self.xs.get()),
                 'y': float(self.ys.get()),
                 'z': float(self.zs.get())}
        try: rx = float(self.xr.get())
        except: rx = self.xr.get()
        try: ry = float(self.yr.get())
        except: ry = self.yr.get()
        try: rz = float(self.zr.get())
        except: rz = self.zr.get()
        rotate = {'x': rx,
                  'y': ry,
                  'z': rz}
        
        # Read current py script and re-write with new settings
        print('Updating settings in: {}, from current GUI settings'.format(self.filename))
        lines = psm.read(self.filename)
        with open(self.filename, 'w') as f:
            inputsflag = True
            for line in lines:
                # if line.startswith('use_GUI') and inputsflag:
                #     line = psm.parse_and_modify(line, True, stringflag=False, splitchar='=')
                # if line.startswith('topofile') and inputsflag:
                #     line = psm.parse_and_modify(line, topofile, stringflag=True, splitchar='=')
                # if line.startswith('nta_file') and inputsflag:
                #     line = psm.parse_and_modify(line, nta_file, stringflag=True, splitchar='=')
                if line.startswith('add2box') and inputsflag:
                    line = psm.parse_and_modify(line, add2box, stringflag=False, splitchar='=')
                if line.startswith('shift') and inputsflag:
                    line = psm.parse_and_modify(line, str(shift), stringflag=False, splitchar='=')
                if line.startswith('rotate') and inputsflag:
                    line = psm.parse_and_modify(line, str(rotate), stringflag=False, splitchar='=')
                if line.startswith('ignore_missing_parameters') and inputsflag:
                    line = psm.parse_and_modify(line, ignore_missing_parameters, stringflag=False, splitchar='=')
                if line.startswith('frc_file') and inputsflag:
                    line = psm.parse_and_modify(line, frc_file, stringflag=True, splitchar='=')
                if line.startswith('assumed') and inputsflag:
                    line = psm.parse_and_modify(line, assumed, stringflag=True, splitchar='=')
                if line.startswith('atom_style') and inputsflag:
                    line = psm.parse_and_modify(line, atom_style, stringflag=True, splitchar='=')
                if line.startswith('parent_directory') and inputsflag:
                    line = psm.parse_and_modify(line, parent_directory, stringflag=True, splitchar='=')
                if line.startswith('newfile') and inputsflag:
                    line = psm.parse_and_modify(line, newfile, stringflag=True, splitchar='=')
                if line.startswith('ff_class') and inputsflag:
                    line = psm.parse_and_modify(line, ff_class, stringflag=True, splitchar='=')
                if line.startswith('use_auto_equivalence') and inputsflag:
                    line = psm.parse_and_modify(line, use_auto_equivalence, stringflag=False, splitchar='=')
                if line.startswith('use_morse_bonds') and inputsflag:
                    line = psm.parse_and_modify(line, use_morse_bonds, stringflag=False, splitchar='=')
                if line.startswith('use_assumed_auto_fill') and inputsflag:
                    line = psm.parse_and_modify(line, use_assumed_auto_fill, stringflag=False, splitchar='=')
                if line.startswith('reset_molids') and inputsflag:
                    line = psm.parse_and_modify(line, reset_molids, stringflag=False, splitchar='=')
                if line.startswith('reset_charges') and inputsflag:
                    line = psm.parse_and_modify(line, reset_charges, stringflag=False, splitchar='=')
                if line.startswith('write_txt_comments') and inputsflag:
                    line = psm.parse_and_modify(line, write_txt_comments, stringflag=False, splitchar='=')
                if line.startswith('write_bond_react') and inputsflag:
                    line = psm.parse_and_modify(line, write_bond_react, stringflag=False, splitchar='=')
                if line.startswith('include_type_labels') and inputsflag:
                    line = psm.parse_and_modify(line, include_type_labels, stringflag=False, splitchar='=')
                if line.startswith('if __name__ == "__main__":'): inputsflag = False
                f.write(line)
        return