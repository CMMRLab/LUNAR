# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
May 13th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.io_functions as io_functions
from src.cell_builder.main import main
import src.py_script_modifier as psm
from tkinter.scrolledtext import ScrolledText
from tkinter import filedialog
from tkinter import Toplevel
from tkinter import ttk
import tkinter as tk
import traceback
import math
import os


############################################
# Function to try to get qty from filepath #
############################################
def get_file_qty(file, delimiter='qty='):
    try:
        # Assume file naming convention is:
        #    *qty=N.data
        # and get tag location in filename
        basename = file[:file.rfind('.')]; qty = 0;
        string_split = basename.split(delimiter)
        qty_guess = string_split[-1].strip()
        
        # Check if qty_guess is valid, if so update qty int
        try:
            qty = int(qty_guess);
            print(f"{file}\nhad *qty=N.data format, where the qty was found to be {qty}\n")
        except: pass
    except: 
        if file.endswith('data'):
            qty = 1
        else: qty = 0
    return qty


##############################
# LUNAR/bond_react_merge GUI #
##############################
class cell_builder_GUI:
    def __init__(self, files, force_field_joining, duplicate, distance_scale, newfile, atom_style, parent_directory, max_rotations,
                 reset_molids, unwrap_atoms_via_image_flags, include_type_labels, group_monomers_locally, seed, domain, GUI_zoom,
                 nfiles=6, scroll_bar=False):
        
        # Find present working directory
        self.pwd = os.getcwd()
        self.filepath = self.pwd
        self.filename = os.path.join(self.pwd, 'cell_builder.py')

        # Initialize window
        self.root = tk.Tk()
        self.root.title('LUNAR/cell_builder.py GUI v1.0')
        #self.root.geometry('600x400')
        
        # Initalize main frame
        if not scroll_bar:
            self.root.resizable(width=False, height=False)
            self.frame = tk.Frame(self.root)
            self.frame.pack()
        
        # Initialize window with a scroll bar
        else:
            GUI_SF = GUI_zoom/100
            height = 25*nfiles + 400
            width = 1250
            height = int(math.ceil(GUI_SF*height))
            width = int(math.ceil(GUI_SF*width))
            if GUI_SF > 1.0: width += int(math.ceil(width/3.25*GUI_SF))
            self.root.minsize(width, height)
            self.frame1 = tk.Frame(self.root)
            self.frame1.pack(fill=tk.BOTH, expand=1)
            self.canvas = tk.Canvas(self.frame1)
            self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
            self.scrollbar = ttk.Scrollbar(self.frame1, orient=tk.VERTICAL, command=self.canvas.yview)
            self.scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
            self.canvas.configure(yscrollcommand=self.scrollbar.set)
            self.canvas.bind('<Configure>', lambda e: self.canvas.configure(scrollregion = self.canvas.bbox('all')))
            self.frame = tk.Frame(self.canvas)
            self.canvas.create_window((0,0), window=self.frame, anchor='nw')
        
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
        self.GUI_zoom = GUI_zoom
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
        self.inputs_frame.grid(row=0, column=0, columnspan=2, padx=xpadding, pady=int(ypadding/2))
        
        # file selection labels
        self.file_label = tk.Label(self.inputs_frame, text='files stack', font=font_settings)
        self.file_label.grid(column=0, row=0, columnspan=2)
        self.file_label = tk.Label(self.inputs_frame, text='qty', font=font_settings)
        self.file_label.grid(column=2, row=0)
        
        # file selection button and qty
        lst_files = list(files.items()); self.nfiles = nfiles; self.files = []; self.qtys = [];
        for n in range(1, self.nfiles+1):
            try: intialfile = list(lst_files)[n-1][0]; intialqty = list(lst_files)[n-1][1];
            except: intialfile = ''; intialqty = '';
            self.file = tk.Entry(self.inputs_frame, width=int(1.35*maxwidth), font=font_settings)
            self.file.insert(0, intialfile)
            self.file.grid(column=0, row=n, columnspan=2)
            self.files.append(self.file)
            self.qty = tk.Entry(self.inputs_frame, width=int(maxwidth/8), font=font_settings)
            self.qty.insert(0, intialqty)
            self.qty.grid(column=2, row=n)
            self.qtys.append(self.qty)
            
        # Button to add a file
        self.file_button = tk.Button(self.inputs_frame, text='add file(s) to stack', font=font_settings, command=self.infile_path)
        self.file_button.grid(column=0, row=self.nfiles+1, columnspan=1)
        
        # Button to remove a file
        self.remove_button = tk.Button(self.inputs_frame, text='remove last file from stack', font=font_settings, width=int(maxwidth/1.85), command=self.remove_file)
        self.remove_button.grid(column=1, row=self.nfiles+1, sticky='news', columnspan=1)
        
        # Button to clear all files
        self.clear_button = tk.Button(self.inputs_frame, text='clear stack', font=font_settings, width=int(maxwidth/13.75), command=self.clear_all)
        self.clear_button.grid(column=2, row=self.nfiles+1, columnspan=1)
        
        # parent_directory entry
        self.parent_directory = tk.Entry(self.inputs_frame, width=int(1.25*maxwidth), font=font_settings)
        self.parent_directory.insert(0, parent_directory)
        self.parent_directory.grid(column=1, row=self.nfiles+2, columnspan=2)
        self.dir_button = tk.Button(self.inputs_frame, text='parent_directory', font=font_settings, command=self.directory_path)
        self.dir_button.grid(column=0, row=self.nfiles+2)
        
        # newfile entry
        self.newfile = tk.Entry(self.inputs_frame, width=int(1.25*maxwidth), font=font_settings)
        self.newfile.insert(0, newfile)
        self.newfile.grid(column=1, row=self.nfiles+3, columnspan=2)
        self.newfile_label = tk.Label(self.inputs_frame, text='newfile', font=font_settings)
        self.newfile_label.grid(column=0, row=self.nfiles+3)

        # Add padding to all frames in self.inputs_frame
        for widget in self.inputs_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))
            
        
        #---------------#
        # Options frame #
        #---------------#
        # Initalize  options frame
        self.options_frame = tk.LabelFrame(self.frame, text='Options', font=font_settings)
        self.options_frame.grid(row=1, column=0, columnspan=2, sticky='news', padx=xpadding, pady=ypadding)
        
        # duplicate entry
        self.duplicate = ttk.Entry(self.options_frame, width=int(maxwidth/8), font=font_settings)
        self.duplicate.insert(0, duplicate)
        self.duplicate.grid(column=0, row=1)
        self.duplicate_label = tk.Label(self.options_frame, text='duplicate', font=font_settings)
        self.duplicate_label.grid(column=0, row=0)
        
        # domain entry
        self.domain = ttk.Entry(self.options_frame, width=int(maxwidth/6), font=font_settings)
        self.domain.insert(0, domain)
        self.domain.grid(column=1, row=1)
        self.domain_label = tk.Label(self.options_frame, text='domain (cubic or Ni x Nj x Nk)', font=font_settings)
        self.domain_label.grid(column=1, row=0)
        
        # distance_scale drop down menu
        self.distance_scale = ttk.Entry(self.options_frame, width=int(maxwidth/11), font=font_settings)
        self.distance_scale.insert(0, distance_scale)
        self.distance_scale.grid(column=2, row=1)
        self.distance_scale_label = tk.Label(self.options_frame, text='distance_scale', font=font_settings)
        self.distance_scale_label.grid(column=2, row=0)
        
        # atom_style drop down
        styles = ['full', 'charge', 'molecular', 'angle', 'bond', 'atomic', 'dipole', 'dpd', 'line']
        self.atom_style = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/11), font=font_settings)
        self.atom_style.current(styles.index(atom_style))
        self.atom_style.grid(column=3, row=1)
        self.atom_style_label = tk.Label(self.options_frame, text='atom_style', font=font_settings)
        self.atom_style_label.grid(column=3, row=0)

        # include_type_labels drop down menu
        styles = [True, False]
        self.include_type_labels = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/11), font=font_settings)
        self.include_type_labels.current(styles.index(include_type_labels))
        self.include_type_labels.grid(column=4, row=1)
        self.include_type_labels_label = tk.Label(self.options_frame, text='include_type_labels', font=font_settings)
        self.include_type_labels_label.grid(column=4, row=0)
        
        # reset_molids drop down menu
        styles = ['files', 'offset', 'clusters', 'skip']
        self.reset_molids = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/11), font=font_settings)
        self.reset_molids.current(styles.index(reset_molids))
        self.reset_molids.grid(column=5, row=1)
        self.reset_molids_label = tk.Label(self.options_frame, text='reset_molids', font=font_settings)
        self.reset_molids_label.grid(column=5, row=0)
        
        # reset_molids drop down menu
        styles = [True, False]
        self.unwrap_atoms_via_image_flags = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/11), font=font_settings)
        self.unwrap_atoms_via_image_flags.current(styles.index(unwrap_atoms_via_image_flags))
        self.unwrap_atoms_via_image_flags.grid(column=6, row=1)
        self.unwrap_atoms_via_image_flags_label = tk.Label(self.options_frame, text='unwrap_atoms_via_image_flags', font=font_settings)
        self.unwrap_atoms_via_image_flags_label.grid(column=6, row=0)
        
        # group_monomers_locally drop down menu
        styles = [True, False]
        self.group_monomers_locally = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/11), font=font_settings)
        self.group_monomers_locally.current(styles.index(group_monomers_locally))
        self.group_monomers_locally.grid(column=0, row=3)
        self.group_monomers_locally_label = tk.Label(self.options_frame, text='group_monomers_locally', font=font_settings)
        self.group_monomers_locally_label.grid(column=0, row=2)
        
        # seed
        self.seed = ttk.Entry(self.options_frame, width=int(maxwidth/6), font=font_settings)
        self.seed.insert(0, seed)
        self.seed.grid(column=1, row=3)
        self.seed_label = tk.Label(self.options_frame, text='seed', font=font_settings)
        self.seed_label.grid(column=1, row=2)
        
        # max X-rotation
        self.mxr = ttk.Entry(self.options_frame, width=int(maxwidth/11), font=font_settings)
        self.mxr.insert(0, max_rotations['x'])
        self.mxr.grid(column=2, row=3)
        self.mxr_label = tk.Label(self.options_frame, text='max X-rotation', font=font_settings)
        self.mxr_label.grid(column=2, row=2)
        
        # max Y-rotation
        self.myr = ttk.Entry(self.options_frame, width=int(maxwidth/11), font=font_settings)
        self.myr.insert(0, max_rotations['y'])
        self.myr.grid(column=3, row=3)
        self.myr_label = tk.Label(self.options_frame, text='max Y-rotation', font=font_settings)
        self.myr_label.grid(column=3, row=2)
        
        # max Z-rotation
        self.mzr = ttk.Entry(self.options_frame, width=int(maxwidth/11), font=font_settings)
        self.mzr.insert(0, max_rotations['z'])
        self.mzr.grid(column=4, row=3)
        self.mzr_label = tk.Label(self.options_frame, text='max Z-rotation', font=font_settings)
        self.mzr_label.grid(column=4, row=2)
        
        # offset_coeff_types drop down menu
        styles = ['none', 'merge', 'offset']
        self.force_field_joining = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/11), font=font_settings)
        self.force_field_joining.current(styles.index(force_field_joining))
        self.force_field_joining.grid(column=5, row=3)
        self.force_field_joining_label = tk.Label(self.options_frame, text='force_field_joining', font=font_settings)
        self.force_field_joining_label.grid(column=5, row=2)


        # Add padding to all frames in self.inputs_frame
        for widget in self.options_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))
            
        
        #------------#
        # Run button #
        #------------#
        self.run = tk.Button(self.frame, text='Run LUNAR/cell_builder.py', font=font_settings, command=self.run_LUNAR)
        self.run.grid(row=2, column=0, columnspan=2, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
                
        #-----------------#
        # update defaults #
        #-----------------#
        self.update = tk.Button(self.frame, text='Save the current GUI settings as the default GUI settings', font=font_settings, command=self.update_py_script)
        self.update.grid(row=3, column=0, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
        
        #------------#
        # Quick help #
        #------------#
        self.quick_help = tk.Button(self.frame, text='Quick help', font=font_settings, command=self.quickhelp)
        self.quick_help.grid(row=3, column=1, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
        
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
            txt = os.path.join(self.pwd, 'src/GUI_quick_help_pages/cell_builder.txt')
            logged = []
            with open(txt, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    # Strip comment's and split by whitespace
                    line = line.split('#')[0]
                    line = line.rstrip()
                    logged.append(line)
        except: # except something failed
            logged.append('FAILED to read LUNAR/src/GUI_quick_help_pages/cell_builder.txt document.')
            logged.append('Most likely cause is the cell_builder.txt file was renamed in the')
            logged.append('LUNAR/src/GUI_quick_help_pages/ directory or directory names were changed.')
        self.popup(logged, title='Quick help')
        return  
    
    # Function to get filepath for topofile
    def infile_path(self):
        ftypes = (('all files', '*.*'), ('data files', '*.data *.data.gz'))
        paths = filedialog.askopenfilenames(initialdir=self.filepath, title='Open file(s)?', filetypes=ftypes)
        if paths:
            for path in paths:
                self.filepath = os.path.dirname(os.path.abspath(path))
                path = os.path.relpath(path)
                
                # Try getting tag from file format naming convention
                qty = get_file_qty(path, delimiter='qty=')
                
                # Try adding file and possibly qty to stack
                try:
                    blanks = [n for n, i in enumerate(self.files) if i.get() == '']
                    self.files[blanks[0]].delete(0, tk.END)
                    self.files[blanks[0]].insert(0, path)
                    if qty > 0: # shortcut naming
                        self.qtys[blanks[0]].delete(0, tk.END)
                        self.qtys[blanks[0]].insert(0, str(qty))
                    elif path.endswith('data'): # avoid zero for grouping file
                        self.qtys[blanks[0]].delete(0, tk.END)
                        self.qtys[blanks[0]].insert(0, '1')
                    else: # final attempt
                        self.qtys[blanks[0]].delete(0, tk.END)
                        self.qtys[blanks[0]].insert(0, '0')
                except: 
                    print('GUI file limit reached. Attempting to add overloaded file.')
                    try:
                        qty = 0
                        if path.endswith('data'):
                            qty = 1
                        self.overloadfile = path; self.overloadqty = qty;
                        self.add_overloaded_filebox()
                    except: print(f'GUI file limit reached and overload add FAILED. Update maxfiles variable in {os.path.relpath(self.filename)}')
        return
    
    # Function to add files to GUI, during overload conditions
    def add_overloaded_filebox(self):
        # adjust based on GUI_SF
        GUI_SF = self.GUI_zoom/100
        font_size = int(math.ceil(GUI_SF*self.font_size))
        xpadding = int(math.ceil(GUI_SF*self.xpadding))
        ypadding = int(math.ceil(GUI_SF*self.ypadding))
        maxwidth = int(math.ceil(GUI_SF*self.maxwidth))
        font_settings = (self.font_type, font_size)
        
        # Add file box
        self.nfiles += 1
        self.file = tk.Entry(self.inputs_frame, width=int(1.35*maxwidth), font=font_settings)
        self.file.insert(0, self.overloadfile)
        self.file.grid(column=0, row=self.nfiles, columnspan=2)
        self.files.append(self.file)
        self.qty = tk.Entry(self.inputs_frame, width=int(maxwidth/8), font=font_settings)
        self.qty.insert(0, self.overloadqty)
        self.qty.grid(column=2, row=self.nfiles)
        self.qtys.append(self.qty)
        
        # adjust packing of other things in inputs frame
        self.file_button.grid(column=0, row=self.nfiles+1, columnspan=1)
        self.remove_button.grid(column=1, row=self.nfiles+1, sticky='news', columnspan=1)
        self.clear_button.grid(column=2, row=self.nfiles+1, columnspan=1)
        self.parent_directory.grid(column=1, row=self.nfiles+2, columnspan=2)
        self.dir_button.grid(column=0, row=self.nfiles+2)
        self.newfile.grid(column=1, row=self.nfiles+3, columnspan=2)
        self.newfile_label.grid(column=0, row=self.nfiles+3)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.inputs_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))
        return
    
    # Function to remove a file from self.files
    def remove_file(self):
        used = [n for n, i in enumerate(self.files) if i.get() != '']
        if used:
            last_add = max(used)
            self.files[last_add].delete(0, tk.END)
            try: self.qtys[last_add].delete(0, tk.END)
            except: pass
        else: print('No files or qtys left to remove')
        return
    
    # Function to clear all files
    def clear_all(self):
        for n, i in enumerate(self.files):
            try: self.files[n].delete(0, tk.END)
            except: pass
        for n, i in enumerate(self.qtys):
            try: self.qtys[n].delete(0, tk.END)
            except: pass
        return
    
    # Function to get directory
    def directory_path(self):
        path =filedialog.askdirectory(initialdir=self.pwd)
        if path:
            path = os.path.relpath(path)
            self.parent_directory.delete(0, tk.END); self.parent_directory.insert(0, path);
        return
    
    # Closing command    
    def closing(self):
        print('Terminating cell_builder GUI'); self.root.destroy();
        return
    
    # Function to run imported LUNAR functions
    def run_LUNAR(self):  
        # Set up log
        log = io_functions.LUNAR_logger()
        valid_inputs = True
        
        # Get information from GUI
        boolean = {'False':False, 'True':True}        
        parent_directory = self.parent_directory.get() 
        atom_style = self.atom_style.get()   
        include_type_labels = boolean[self.include_type_labels.get()]   
        reset_molids = self.reset_molids.get()
        unwrap_atoms_via_image_flags = boolean[self.unwrap_atoms_via_image_flags.get()]
        try: duplicate = int(self.duplicate.get())
        except:
            log.GUI_error('ERROR duplicate is not an int')
            valid_inputs = False
        try: distance_scale = float(self.distance_scale.get())
        except:
            log.GUI_error('ERROR distance_scale is not a float')
            valid_inputs = False
        try:
            max_rotations = {'x': float(self.mxr.get()),
                             'y': float(self.myr.get()),
                             'z': float(self.mzr.get())}
        except:
            log.GUI_error('ERROR max X, Y, or Z rotation is not a float')
            valid_inputs = False
        try: seed = int(self.seed.get())
        except:
            log.GUI_error('ERROR seed is not an int')
            valid_inputs = False
        newfile = self.newfile.get()
        group_monomers_locally = boolean[self.group_monomers_locally.get()]
        force_field_joining = self.force_field_joining.get()
        domain = self.domain.get()
        
        # build files based on GUI inputs
        files = {} # {file/path : qty }
        for i, j  in zip(self.files, self.qtys):
            try: 
                if j.get() != '' and i.get() != '': files[i.get()] = int(j.get())
                elif j.get() != '' or i.get() != '': log.warn(f'WARNING incomplete information specified from GUI. file: {i.get()} qty: {j.get()}') 
            except: pass

        # Run LUNAR/cell_builder on a different thread
        if valid_inputs:
            try: main(files, force_field_joining, duplicate, distance_scale, newfile, atom_style, parent_directory, max_rotations,
                      reset_molids, unwrap_atoms_via_image_flags, include_type_labels, group_monomers_locally, seed, domain, log=log)
            except Exception:
                log.GUI_error(traceback.format_exc())
        self.popup(log.logged, title='Outputs', width=200)
        return   
    
    # Function to pop-up scrollable text
    def popup(self, out, title='Outputs', width=150):
        page = Toplevel(self.root)
        page.title(title)
        outputs = ScrolledText(page, height=30, width=width, font=('consolas', '12', 'normal'))
        outputs.pack()
        outputs.insert(tk.INSERT, '\n'.join(out))
        outputs.config(state=tk.DISABLED)
        page.resizable(width=False, height=False)
        return
    
    # Function to update py script default settings
    def update_py_script(self):
        # Get information from GUI
        boolean = {'False':False, 'True':True}        
        parent_directory = io_functions.path_to_string(self.parent_directory.get()) 
        atom_style = self.atom_style.get()   
        include_type_labels = boolean[self.include_type_labels.get()]   
        reset_molids = self.reset_molids.get()
        unwrap_atoms_via_image_flags = boolean[self.unwrap_atoms_via_image_flags.get()]
        duplicate = int(self.duplicate.get())
        distance_scale = float(self.distance_scale.get())
        newfile = self.newfile.get()
        max_rotations = {'x': float(self.mxr.get()),
                         'y': float(self.myr.get()),
                         'z': float(self.mzr.get())}
        group_monomers_locally = boolean[self.group_monomers_locally.get()]
        seed = int(self.seed.get())
        force_field_joining = self.force_field_joining.get()
        domain = self.domain.get()
        
        # Read current py script and re-write with new settings
        print('Updating settings in: {}, from current GUI settings'.format(self.filename))
        lines = psm.read(self.filename)
        with open(self.filename, 'w') as f:
            inputsflag = True
            for line in lines:
                # if line.startswith('use_GUI') and inputsflag:
                #     line = psm.parse_and_modify(line, True, stringflag=False, splitchar='=')
                if line.startswith('parent_directory') and inputsflag:
                    line = psm.parse_and_modify(line, parent_directory, stringflag=True, splitchar='=')
                if line.startswith('domain') and inputsflag:
                    line = psm.parse_and_modify(line, domain, stringflag=True, splitchar='=')
                if line.startswith('group_monomers_locally') and inputsflag:
                    line = psm.parse_and_modify(line, group_monomers_locally, stringflag=False, splitchar='=')
                if line.startswith('atom_style') and inputsflag:
                    line = psm.parse_and_modify(line, atom_style, stringflag=True, splitchar='=')
                if line.startswith('include_type_labels') and inputsflag:
                    line = psm.parse_and_modify(line, include_type_labels, stringflag=False, splitchar='=')
                if line.startswith('reset_molids') and inputsflag:
                    line = psm.parse_and_modify(line, reset_molids, stringflag=True, splitchar='=')
                if line.startswith('unwrap_atoms_via_image_flags') and inputsflag:
                    line = psm.parse_and_modify(line, unwrap_atoms_via_image_flags, stringflag=False, splitchar='=')
                if line.startswith('force_field_joining') and inputsflag:
                    line = psm.parse_and_modify(line, force_field_joining, stringflag=False, splitchar='=')
                if line.startswith('duplicate') and inputsflag:
                    line = psm.parse_and_modify(line, duplicate, stringflag=False, splitchar='=')
                if line.startswith('seed') and inputsflag:
                    line = psm.parse_and_modify(line, seed, stringflag=False, splitchar='=')
                if line.startswith('distance_scale') and inputsflag:
                    line = psm.parse_and_modify(line, distance_scale, stringflag=False, splitchar='=')
                if line.startswith('newfile') and inputsflag:
                    line = psm.parse_and_modify(line, newfile, stringflag=True, splitchar='=')
                if line.startswith('max_rotations') and inputsflag:
                    line = psm.parse_and_modify(line, str(max_rotations), stringflag=False, splitchar='=')
                if line.startswith('if __name__ == "__main__":'): inputsflag = False
                f.write(line)
        return