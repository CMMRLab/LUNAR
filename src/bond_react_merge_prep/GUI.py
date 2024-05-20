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
from src.bond_react_merge_prep.main import main
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


###################################
# LUNAR/bond_react_merge_prep GUI #
###################################
class bond_react_merge_prep_GUI:
    def __init__(self, topofile, cta_file, newfile, atom_style, parent_directory, rm_unused_coeffs, GUI_zoom):
        
        # Find present working directory
        self.pwd = os.getcwd()
        self.filepath = self.pwd
        self.filename = os.path.join(self.pwd, 'bond_react_merge_prep.py')

        # Initialize window
        self.root = tk.Tk()
        self.root.title('LUNAR/bond_react_merge_prep.py GUI v1.0')
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
        self.maxwidth = 100
        
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
        self.inputs_frame.grid(row=0, column=0, padx=xpadding, pady=ypadding)
        
        # topofile selection button
        self.topofile = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.topofile.insert(0, topofile)
        self.topofile.grid(column=1, row=0)
        self.topofile_button = tk.Button(self.inputs_frame, text='topofile', font=font_settings, command=self.topofile_path)
        self.topofile_button.grid(column=0, row=0)
        
        # nta_file selection button
        self.cta_file = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.cta_file.insert(0, cta_file)
        self.cta_file.grid(column=1, row=1)
        self.cta_file_button = tk.Button(self.inputs_frame, text='cta_file', font=font_settings, command=self.cta_file_path)
        self.cta_file_button.grid(column=0, row=1)
        
        # parent_directory entry
        self.parent_directory = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.parent_directory.insert(0, parent_directory)
        self.parent_directory.grid(column=1, row=2)
        self.dir_button = tk.Button(self.inputs_frame, text='parent_directory', font=font_settings, command=self.directory_path)
        self.dir_button.grid(column=0, row=2)
        
        # newfile entry
        self.newfile = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.newfile.insert(0, newfile)
        self.newfile.grid(column=1, row=5)
        self.newfile_label = tk.Label(self.inputs_frame, text='newfile', font=font_settings)
        self.newfile_label.grid(column=0, row=5)
        
        # atom_style drop down
        styles = ['full', 'charge', 'molecular', 'angle', 'bond', 'atomic', 'dipole', 'dpd', 'line']
        self.atom_style = ttk.Combobox(self.inputs_frame, values=styles, width=int(maxwidth/1.02), font=font_settings)
        self.atom_style.current(styles.index(atom_style))
        self.atom_style.grid(column=1, row=6)
        self.atom_style_label = tk.Label(self.inputs_frame, text='atom_style', font=font_settings)
        self.atom_style_label.grid(column=0, row=6)
        
        # include_rcut drop down menu
        styles = [True, False]
        self.rm_unused_coeffs = ttk.Combobox(self.inputs_frame, values=styles, width=int(maxwidth/1.02), font=font_settings)
        self.rm_unused_coeffs.current(styles.index(rm_unused_coeffs))
        self.rm_unused_coeffs.grid(column=1, row=7)
        self.rm_unused_coeffs_label = tk.Label(self.inputs_frame, text='rm_unused_coeffs', font=font_settings)
        self.rm_unused_coeffs_label.grid(column=0, row=7)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.inputs_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=ypadding)
            
        
        #------------#
        # Run button #
        #------------#
        self.run = tk.Button(self.frame, text='Run LUNAR/bond_react_merge_prep.py', font=font_settings, command=self.run_LUNAR)
        self.run.grid(row=2, column=0, sticky='news', padx=xpadding, pady=ypadding)
        
        #-----------------#
        # update defaults #
        #-----------------#
        self.update = tk.Button(self.frame, text='Save the current GUI settings as the default GUI settings', font=font_settings, command=self.update_py_script)
        self.update.grid(row=3, column=0, sticky='news', padx=xpadding, pady=ypadding)
        
        
        #------------------------#
        # Run mainloop and close #
        #------------------------#
        self.root.protocol('WM_DELETE_WINDOW', self.closing)
        self.root.mainloop()
        
        
        
    #################################
    # Functions to call as commands #
    #################################
    # Function to get filepath for topofile
    def topofile_path(self):
        ftypes = (('all files', '*.*'), ('data files', '*.data *.data.gz'), ('mol files', '*.mol'), ('mol2 files', '*.mol2'), ('mdf files', '*.mdf'), ('sdf files', '*.sdf'))
        path = filedialog.askopenfilename(initialdir=self.filepath, title='Open topofile?', filetypes=ftypes)
        if path:
            self.filepath = os.path.dirname(os.path.abspath(path))
            path = os.path.relpath(path)
            self.topofile.delete(0, tk.END); self.topofile.insert(0, path);
        return
    
    # Function to get directory
    def directory_path(self):
        path =filedialog.askdirectory(initialdir=self.pwd)
        if path:
            path = os.path.relpath(path)
            self.parent_directory.delete(0, tk.END); self.parent_directory.insert(0, path);
        return
    
    # Function to get filepath for cta_file
    def cta_file_path(self):
        ftypes = (('all files', '*.*'), ('cta files', '*.nta'), ('car files', '*.car'))
        path = filedialog.askopenfilename(initialdir=self.filepath, title='Open nta_file?', filetypes=ftypes)
        if path:
            self.filepath = os.path.dirname(os.path.abspath(path))
            path = os.path.relpath(path)
            self.cta_file.delete(0, tk.END); self.cta_file.insert(0, path);
        return
    
    # Closing command    
    def closing(self):
        print('Terminating bond_react_merge_prep GUI'); self.root.destroy();
        return
    
    # Function to run imported LUNAR functions
    def run_LUNAR(self):  
        # Set up log
        log = io_functions.LUNAR_logger()
        valid_inputs = True
        
        # Get information from GUI
        boolean = {'False':False, 'True':True}
        topofile = self.topofile.get()
        cta_file = self.cta_file.get()
        parent_directory = self.parent_directory.get() 
        atom_style = self.atom_style.get()
        newfile = self.newfile.get()
        rm_unused_coeffs = boolean[self.rm_unused_coeffs.get()]

        # Run LUNAR/bond_react_merge_prep
        if valid_inputs:
            try: 
                inputs = (topofile, cta_file, newfile, atom_style, parent_directory, rm_unused_coeffs, [], log)
                
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
        topofile = io_functions.path_to_string(self.topofile.get())
        cta_file = io_functions.path_to_string(self.cta_file.get())
        parent_directory = io_functions.path_to_string(self.parent_directory.get()) 
        atom_style = self.atom_style.get()
        newfile = self.newfile.get()
        rm_unused_coeffs = boolean[self.rm_unused_coeffs.get()]
        
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
                # if line.startswith('cta_file') and inputsflag:
                #     line = psm.parse_and_modify(line, cta_file, stringflag=True, splitchar='=')
                if line.startswith('atom_style') and inputsflag:
                    line = psm.parse_and_modify(line, atom_style, stringflag=True, splitchar='=')
                if line.startswith('parent_directory') and inputsflag:
                    line = psm.parse_and_modify(line, parent_directory, stringflag=True, splitchar='=')
                if line.startswith('newfile') and inputsflag:
                    line = psm.parse_and_modify(line, newfile, stringflag=True, splitchar='=')
                if line.startswith('rm_unused_coeffs') and inputsflag:
                    line = psm.parse_and_modify(line, rm_unused_coeffs, stringflag=False, splitchar='=')
                if line.startswith('if __name__ == "__main__":'): inputsflag = False
                f.write(line)
        return