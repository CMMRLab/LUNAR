# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.3
May 13th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.GUI_scale_settings as GUI_scale_settings
from src.bond_react_merge.main import main
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


############################################
# Function to try to get tag from filepath #
############################################
def get_file_tag(file, GUI_log, delimiter='tag='):
    try:
        # Assume file naming convention is:
        #    *tag=filetag.data
        # and get tag location in filename
        basename = file[:file.rfind('.')]; tag = '';
        string_split = basename.split(delimiter)
        tag_guess = string_split[-1].strip()
        
        # strip last number of tag_guess. Assumes tags either:
        # - data1, data2, ... dataN  ->  data
        # - pre1, pre2, ... preN     ->  pre
        # - post1, post2, ... postN  ->  post
        # - infile                   ->  infile
        tmpname = ''.join([i for i in tag_guess if i.isalpha()])
        
        # Check if tmpname is valid, if so update tag string
        if tmpname in ['pre', 'post', 'data', 'infile']:
            tag = tag_guess
            GUI_log.out(f"{file}\nhad *tag=filetag.data format, where the tag was found to be {tag}\n")
    except: tag = '';
    return tag

##########################################################
# Function to guess file tag from info in filename (most #
# users will be naming files in pre/post with an int)    #
##########################################################
def guess_file_tag(file, ntags, GUI_log):
    root = os.path.basename(file)
    #root = root[:root.rfind('.')]
    if 'pre' in root or 'post' in root or 'data' in root:
        if 'pre' in root: tag_guess = 'pre'
        elif 'post' in root: tag_guess = 'post'
        elif 'data' in root: tag_guess = 'data'
        else: tag_guess = 'data'
        int_strings = [[]]; count = 0;
        for n, i in enumerate(root):
            try: next_char = root[n+1]
            except: next_char = ''
            if i.isdigit():
                if not next_char.isdigit():
                    int_strings[count].append(i)
                else:
                    count += 1
                    int_strings.append([])
                    int_strings[count] = [i]
        try: 
            ints = [int(''.join(i)) for i in int_strings if i]
            ID = str(min(ints))
            ntags[tag_guess] += 1
            tag_guess = tag_guess + ID
        except:
            ID = str(ntags[tag_guess]+1)
            ntags[tag_guess] += 1
            tag_guess = tag_guess + ID
    else:
        if 'txt' in root or 'merge' in root:
            tag_guess = 'infile'
        else:
            tag_guess = 'data' + str(ntags['data']+1)
            ntags['data'] += 1  
    GUI_log.out(f"{file}\n was guessed to have tag {tag_guess}. PLEASE REVIEW GUESSED TAG!\n")
    return tag_guess


##############################
# LUNAR/bond_react_merge GUI #
##############################
class bond_react_merge_GUI:
    def __init__(self, files, parent_directory, newfile, atom_style, generate_map_file, write_rxn_mol2files, 
                 write_rxn_datafiles, write_moleculefiles, print_options, commandline_inputs, map_near_edge_rxn_charges,
                 molecule_file_options, include_type_labels, GUI_zoom, nfiles=6, scroll_bar=False):

        
        # Find present working directory
        self.pwd = os.getcwd()
        self.filename = os.path.join(self.pwd, 'bond_react_merge.py')

        # Initialize window
        self.root = tk.Tk()
        self.root.title('LUNAR/bond_react_merge.py GUI v1.0')
        #self.root.geometry('1000x600') # widthxheight
        
        # Setup GUI log
        self.GUI_log = io_functions.LUNAR_logger()
        self.GUI_log.configure(level='production')

        
        # Initalize main frame
        if not scroll_bar:
            self.root.resizable(width=False, height=False)
            self.frame = tk.Frame(self.root)
            self.frame.pack()

        
        # Initialize window with a scroll bar
        else:
            GUI_SF = GUI_zoom/100
            height = 20*nfiles + 450
            width = 1000
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
        self.maxwidth = 125
        
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
        self.maxwidth = 125

        # Check if user specified any other nong-global scaling settings
        #   extra maxwidth scaling becasuse this is a wide GUI
        scale_settings = GUI_scale_settings.screen_settings
        if 'scaling_factor' in scale_settings:
            if isinstance(scale_settings['scaling_factor'], (int, float)):
               self.xpadding = int(self.xpadding/scale_settings['scaling_factor'])
               self.ypadding = int(self.ypadding/scale_settings['scaling_factor'])
               self.maxwidth = int(self.maxwidth/(scale_settings['scaling_factor']*1.2))
        
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
        self.inputs_frame.grid(row=0, column=0, columnspan=2, padx=xpadding, pady=ypadding)
        
        # file selection labels
        self.file_label = tk.Label(self.inputs_frame, text='files stack', font=font_settings)
        self.file_label.grid(column=0, row=0, columnspan=2)
        self.file_label = tk.Label(self.inputs_frame, text='file-tag', font=font_settings)
        self.file_label.grid(column=2, row=0)
        
        # file selection button and qty
        lst_files = list(files.items()); self.nfiles = nfiles; self.files = []; self.tags = []; self.ntabs = 0;
        for n in range(1, self.nfiles+1):
            try: intialfile = list(lst_files)[n-1][1]; intialtag = list(lst_files)[n-1][0];
            except: intialfile = ''; intialtag = '';
            self.file = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
            self.file.insert(0, intialfile)
            self.file.grid(column=0, row=n, columnspan=2)
            self.files.append(self.file)
            self.tag = tk.Entry(self.inputs_frame, width=int(maxwidth/11.5), font=font_settings)
            self.tag.insert(0, intialtag)
            self.tag.grid(column=2, row=n)
            self.tags.append(self.tag)
        
        # Button to add a file
        self.file_button = tk.Button(self.inputs_frame, text='add file(s) to stack', font=font_settings, command=self.infile_path)
        self.file_button.grid(column=0, row=self.nfiles+1, columnspan=1)
        
        # Button to remove a file
        self.remove_button = tk.Button(self.inputs_frame, text='remove last file from stack', font=font_settings, width=int(maxwidth/1.915), command=self.remove_file)
        self.remove_button.grid(column=1, row=self.nfiles+1, sticky='news', columnspan=1)
        
        # Button to clear all files
        self.clear_button = tk.Button(self.inputs_frame, text='clear stack', font=font_settings, width=int(maxwidth/12), command=self.clear_all)
        self.clear_button.grid(column=2, row=self.nfiles+1, columnspan=1)
        
        # parent_directory entry
        self.parent_directory = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.parent_directory.insert(0, parent_directory)
        self.parent_directory.grid(column=1, row=self.nfiles+2, columnspan=2)
        self.dir_button = tk.Button(self.inputs_frame, text='parent_directory', font=font_settings, command=self.directory_path)
        self.dir_button.grid(column=0, row=self.nfiles+2)
        
        # newfile entry
        self.newfile = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
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
        
        # atom_style drop down
        styles = ['full', 'charge', 'molecular', 'angle', 'bond', 'atomic', 'dipole', 'dpd', 'line']
        self.atom_style = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/8), font=font_settings)
        self.atom_style.current(styles.index(atom_style))
        self.atom_style.grid(column=0, row=1)
        self.atom_style_label = tk.Label(self.options_frame, text='atom_style', font=font_settings)
        self.atom_style_label.grid(column=0, row=0)
        
        # generate_map_file drop down menu
        styles = [True, False]
        self.generate_map_file = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/8), font=font_settings)
        self.generate_map_file.current(styles.index(generate_map_file))
        self.generate_map_file.grid(column=1, row=1)
        self.generate_map_file_label = tk.Label(self.options_frame, text='generate_map_file', font=font_settings)
        self.generate_map_file_label.grid(column=1, row=0)
        
        # write_rxn_mol2files drop down menu
        styles = [True, False]
        self.write_rxn_mol2files = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/8), font=font_settings)
        self.write_rxn_mol2files.current(styles.index(write_rxn_mol2files))
        self.write_rxn_mol2files.grid(column=2, row=1)
        self.write_rxn_mol2files_labels_label = tk.Label(self.options_frame, text='write_rxn_mol2files', font=font_settings)
        self.write_rxn_mol2files_labels_label.grid(column=2, row=0)
        
        # write_rxn_datafiles drop down menu
        styles = [True, False]
        self.write_rxn_datafiles = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/8), font=font_settings)
        self.write_rxn_datafiles.current(styles.index(write_rxn_datafiles))
        self.write_rxn_datafiles.grid(column=3, row=1)
        self.write_rxn_datafiles_label = tk.Label(self.options_frame, text='write_rxn_datafiles', font=font_settings)
        self.write_rxn_datafiles_label.grid(column=3, row=0)
        
        # write_mol_datafiles drop down menu
        styles = [True, False]
        self.write_moleculefiles = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/8), font=font_settings)
        self.write_moleculefiles.current(styles.index(write_moleculefiles))
        self.write_moleculefiles.grid(column=4, row=1)
        self.write_moleculefiles_label = tk.Label(self.options_frame, text='write_moleculefiles', font=font_settings)
        self.write_moleculefiles_label.grid(column=4, row=0)
        
        # map_near_edge_rxn_charges drop down menu
        self.edgeq = [i for i in range(10)]; self.edgeq.insert(0, False);
        self.edgeqmap = {str(i):i for i in self.edgeq}
        self.map_near_edge_rxn_charges = ttk.Combobox(self.options_frame, values=self.edgeq, width=int(maxwidth/10), font=font_settings)
        self.map_near_edge_rxn_charges.current(self.edgeq.index(include_type_labels))
        self.map_near_edge_rxn_charges.grid(column=0, row=3)
        self.map_near_edge_rxn_charges_label = tk.Label(self.options_frame, text='map_near_edge_rxn_charges', font=font_settings)
        self.map_near_edge_rxn_charges_label.grid(column=0, row=2)

        # include_type_labels drop down menu
        styles = [True, False]
        self.include_type_labels = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/10), font=font_settings)
        self.include_type_labels.current(styles.index(include_type_labels))
        self.include_type_labels.grid(column=1, row=3)
        self.include_type_labels_label = tk.Label(self.options_frame, text='include_type_labels', font=font_settings)
        self.include_type_labels_label.grid(column=1, row=2)
                
        # molecule_file_options entry
        self.molecule_file_options = tk.Entry(self.options_frame, width=int(maxwidth/1.5), font=font_settings)
        self.molecule_file_options.insert(0, ','.join(molecule_file_options))
        self.molecule_file_options.grid(column=2, row=3, columnspan=4)
        self.molecule_file_options_label = tk.Label(self.options_frame, text='molecule_file_options (comma separated w/no whitespace)', font=font_settings)
        self.molecule_file_options_label.grid(column=2, row=2, columnspan=4)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.options_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))
            
        
        #------------#
        # Run button #
        #------------#
        self.run = tk.Button(self.frame, text='Run LUNAR/bond_react_merge.py', font=font_settings, command=self.run_LUNAR)
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
            txt = os.path.join(self.pwd, 'src/GUI_quick_help_pages/bond_react_merge.txt')
            logged = []
            with open(txt, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    # Strip comment's and split by whitespace
                    line = line.split('#')[0]
                    line = line.rstrip()
                    logged.append(line)
        except: # except something failed
            logged.append('FAILED to read LUNAR/src/GUI_quick_help_pages/bond_react_merge.txt document.')
            logged.append('Most likely cause is the bond_react_merge.txt file was renamed in the')
            logged.append('LUNAR/src/GUI_quick_help_pages/ directory or directory names were changed.')
        self.popup(logged, title='Quick help')
        return        

    # Function to get filepath for topofile
    def infile_path(self):
        ftypes = (('all files', '*.*'), ('data files', '*.data  *.data.gz'), ('txt files', '*.txt'))
        paths = filedialog.askopenfilenames(title='Open files(s)?', filetypes=ftypes)
        if paths:
            ntags = {'data':0, 'pre':0, 'post':0}
            self.GUI_log.clear_all()
            for path in paths:
                path = os.path.relpath(path)
                
                # Try getting tag from file format naming convention
                tag = get_file_tag(path, self.GUI_log, delimiter='tag=')
                
                # if tag could not be found through naming convention, guess the tag
                if tag == '': tag = guess_file_tag(path, ntags, self.GUI_log)
    
                # Try adding file and possibly tag to stack
                try:
                    blanks = [n for n, i in enumerate(self.files) if i.get() == '']
                    self.files[blanks[0]].delete(0, tk.END)
                    self.files[blanks[0]].insert(0, path)
                    if tag != '':
                        try:
                            self.tags[blanks[0]].delete(0, tk.END)
                            self.tags[blanks[0]].insert(0, tag)
                        except: pass
                except: 
                    self.GUI_log.out('GUI file limit reached. Attempting to add overloaded file.')
                    try:
                        self.overloadfile = path; self.overloadtag = tag;
                        self.add_overloaded_filebox()
                    except: self.GUI_log.out(f'GUI file limit reached and overload add FAILED. Update maxfiles variable in {os.path.relpath(self.filename)}')
            self.popup(self.GUI_log.logged, title='File-tags')
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
        self.file = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.file.insert(0, self.overloadfile)
        self.file.grid(column=0, row=self.nfiles, columnspan=2)
        self.files.append(self.file)
        self.tag = tk.Entry(self.inputs_frame, width=int(maxwidth/11.5), font=font_settings)
        self.tag.insert(0, self.overloadtag)
        self.tag.grid(column=2, row=self.nfiles)
        self.tags.append(self.tag)
        
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
            try: self.tags[last_add].delete(0, tk.END)
            except: pass
        else: print('No files or tags left to remove')
        return
    
    # Function to clear all files
    def clear_all(self):
        for n, i in enumerate(self.files):
            try: self.files[n].delete(0, tk.END)
            except: pass
        for n, i in enumerate(self.tags):
            try: self.tags[n].delete(0, tk.END)
            except: pass
        return
    
    # Function to get directory
    def directory_path(self):
        path =filedialog.askdirectory()
        if path:
            path = os.path.relpath(path)
            self.parent_directory.delete(0, tk.END); self.parent_directory.insert(0, path);
        return

    # Closing command    
    def closing(self):
        print('Terminating bond_react_merge GUI'); self.root.destroy();
        return
    
    # Function to run imported LUNAR functions
    def run_LUNAR(self):  
        # Set up log
        log = io_functions.LUNAR_logger()
        log.configure(level='production', print2console=True,  write2log=True)
        valid_inputs = True
        
        # Get information from GUI
        boolean = {'False':False, 'True':True}
        parent_directory = self.parent_directory.get() 
        newfile = self.newfile.get()
        atom_style = self.atom_style.get()
        generate_map_file = boolean[self.generate_map_file.get()]   
        include_type_labels = boolean[self.include_type_labels.get()]   
        write_rxn_mol2files = boolean[self.write_rxn_mol2files.get()]   
        write_rxn_datafiles = boolean[self.write_rxn_datafiles.get()] 
        write_moleculefiles = boolean[self.write_moleculefiles.get()] 
        try: molecule_file_options = self.molecule_file_options.get().split(',')
        except: molecule_file_options = []
        map_near_edge_rxn_charges = self.edgeqmap[self.map_near_edge_rxn_charges.get()]
        print_options = False
        commandline_inputs = []
        
        # build files based on GUI inputs
        files = {} # { "file-tag" :   'name-of-tagged-datafile' }
        for i, j  in zip(self.files, self.tags):
            try: 
                if j.get() != '' and i.get() != '': files[j.get()] = i.get()
                elif j.get() != '' or i.get() != '': log.out(f'WARNING incomplete information specified from GUI. file-tag: {j.get()} file: {i.get()}') 
            except: pass

        # Run LUNAR/bond_react_merge
        if valid_inputs:
            try:
                inputs = (files, parent_directory, newfile, atom_style, generate_map_file, write_rxn_mol2files, 
                          write_rxn_datafiles, write_moleculefiles, print_options, commandline_inputs, map_near_edge_rxn_charges,
                          molecule_file_options, include_type_labels, log)
                
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
        tk2ff = {'0':0, '1':1, '2':2, 'r':'r', 'd':'d', 's1':'s1', 's2':'s2'}
        boolean = {'False':False, 'True':True}
        parent_directory = io_functions.path_to_string(self.parent_directory.get()) 
        newfile = self.newfile.get()
        atom_style = self.atom_style.get()
        generate_map_file = boolean[self.generate_map_file.get()]   
        include_type_labels = boolean[self.include_type_labels.get()]   
        write_rxn_mol2files = boolean[self.write_rxn_mol2files.get()]   
        write_rxn_datafiles = boolean[self.write_rxn_datafiles.get()] 
        write_moleculefiles = boolean[self.write_moleculefiles.get()] 
        try: molecule_file_options = self.molecule_file_options.get().split(',')
        except: molecule_file_options = []
        map_near_edge_rxn_charges = self.edgeqmap[self.map_near_edge_rxn_charges.get()]
        
        # build files based on GUI inputs
        files = {} # { "file-tag" :   'name-of-tagged-datafile' }
        for i, j  in zip(self.files, self.tags):
            try: 
                if j.get() != '' and i.get() != '': files[j.get()] = i.get()
                elif j.get() != '' or i.get() != '': print(f'WARNING incomplete information specified from GUI. file-tag: {j.get()} file: {i.get()}') 
            except: pass
        
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
                if line.startswith('newfile') and inputsflag:
                    line = psm.parse_and_modify(line, newfile, stringflag=True, splitchar='=') 
                if line.startswith('atom_style') and inputsflag:
                    line = psm.parse_and_modify(line, atom_style, stringflag=True, splitchar='=')
                if line.startswith('generate_map_file') and inputsflag:
                    line = psm.parse_and_modify(line, generate_map_file, stringflag=False, splitchar='=')   
                if line.startswith('include_type_labels') and inputsflag:
                    line = psm.parse_and_modify(line, include_type_labels, stringflag=False, splitchar='=')
                if line.startswith('write_rxn_mol2files') and inputsflag:
                    line = psm.parse_and_modify(line, write_rxn_mol2files, stringflag=False, splitchar='=')
                if line.startswith('write_rxn_datafiles') and inputsflag:
                    line = psm.parse_and_modify(line, write_rxn_datafiles, stringflag=False, splitchar='=')
                if line.startswith('write_moleculefiles') and inputsflag:
                    line = psm.parse_and_modify(line, write_moleculefiles, stringflag=False, splitchar='=')
                if line.startswith('molecule_file_options') and inputsflag:
                    line = psm.parse_and_modify(line, str(molecule_file_options), stringflag=False, splitchar='=')
                if line.startswith('map_near_edge_rxn_charges') and inputsflag:
                    line = psm.parse_and_modify(line, map_near_edge_rxn_charges, stringflag=False, splitchar='=')
                if line.startswith('if __name__ == "__main__":'): inputsflag = False
                f.write(line)
        return