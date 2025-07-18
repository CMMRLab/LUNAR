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
from src.add_pi_electrons.main import main
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


# Function to convert strings2strings and ints2ints
def s2s_i2i(inlst):
    outlst = []
    for i in inlst:
        try: j = int(i)
        except: j = str(i)
        outlst.append(j)
    return outlst


##############################
# LUNAR/add_pi_electrons GUI #
##############################
class add_pi_electrons_GUI:
    def __init__(self, topofile, types2convert, atom_style, reset_charges, net_zero_charge, convert2cg1, add_pi_electrons,
                 parent_directory, newfile, include_type_labels, neighbor_charge_constraint, reset_simulation_cell, GUI_zoom):
        
        # Find present working directory
        self.pwd = os.getcwd()
        self.filename = os.path.join(self.pwd, 'add_pi_electrons.py')

        # Initialize window
        self.root = tk.Tk()
        self.root.title('LUNAR/add_pi_electrons.py GUI v1.0')
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
        self.maxwidth = 85
        
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
        self.maxwidth = 85

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
        self.topofile = tk.Entry(self.inputs_frame, width=int(1.4*maxwidth), font=font_settings)
        self.topofile.insert(0, topofile)
        self.topofile.grid(column=1, row=0)
        self.topofile_button = tk.Button(self.inputs_frame, text='topofile', font=font_settings, command=self.topofile_path)
        self.topofile_button.grid(column=0, row=0)

        # parent_directory entry
        self.parent_directory = tk.Entry(self.inputs_frame, width=int(1.4*maxwidth), font=font_settings)
        self.parent_directory.insert(0, parent_directory)
        self.parent_directory.grid(column=1, row=4)
        self.dir_button = tk.Button(self.inputs_frame, text='parent_directory', font=font_settings, command=self.directory_path)
        self.dir_button.grid(column=0, row=4)
        
        # newfile entry
        self.newfile = tk.Entry(self.inputs_frame, width=int(1.4*maxwidth), font=font_settings)
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
                

        
        # atom_style drop down
        styles = ['full', 'charge', 'molecular', 'angle', 'bond', 'atomic', 'dipole', 'dpd', 'line']
        self.atom_style = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/8.5), font=font_settings)
        self.atom_style.current(styles.index(atom_style))
        self.atom_style.grid(column=0, row=1)
        self.atom_style_label = tk.Label(self.options_frame, text='atom_style', font=font_settings)
        self.atom_style_label.grid(column=0, row=0)
        
        # convert2cg1 drop down menu
        styles = [True, False]
        self.convert2cg1 = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/8.5), font=font_settings)
        self.convert2cg1.current(styles.index(convert2cg1))
        self.convert2cg1.grid(column=1, row=1)
        self.convert2cg1_label = tk.Label(self.options_frame, text='convert2cg1', font=font_settings)
        self.convert2cg1_label.grid(column=1, row=0)
        
        # add_pi_electrons drop down menu
        styles = [True, False]
        self.add_pi_electrons = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/8.5), font=font_settings)
        self.add_pi_electrons.current(styles.index(add_pi_electrons))
        self.add_pi_electrons.grid(column=2, row=1)
        self.add_pi_electrons_label = tk.Label(self.options_frame, text='add_pi_electrons', font=font_settings)
        self.add_pi_electrons_label.grid(column=2, row=0)
        
        # reset_charges drop down menu
        styles = [True, False]
        self.reset_charges = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/8.5), font=font_settings)
        self.reset_charges.current(styles.index(reset_charges))
        self.reset_charges.grid(column=3, row=1)
        self.reset_charges_label = tk.Label(self.options_frame, text='reset_charges', font=font_settings)
        self.reset_charges_label.grid(column=3, row=0)
        
        # net_zero_charge drop down menu
        styles = [True, False]
        self.net_zero_charge = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/8.5), font=font_settings)
        self.net_zero_charge.current(styles.index(net_zero_charge))
        self.net_zero_charge.grid(column=4, row=1)
        self.net_zero_charge_label = tk.Label(self.options_frame, text='net_zero_charge', font=font_settings)
        self.net_zero_charge_label.grid(column=4, row=0)
        
        # include_type_labels drop down menu
        styles = [True, False]
        self.reset_simulation_cell = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/8.5), font=font_settings)
        self.reset_simulation_cell.current(styles.index(reset_simulation_cell))
        self.reset_simulation_cell.grid(column=5, row=1)
        self.reset_simulation_cell_label = tk.Label(self.options_frame, text='reset_simulation_cell', font=font_settings)
        self.reset_simulation_cell_label.grid(column=5, row=0)
        
        # include_type_labels drop down menu
        styles = [True, False]
        self.include_type_labels = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/8.5), font=font_settings)
        self.include_type_labels.current(styles.index(include_type_labels))
        self.include_type_labels.grid(column=0, row=3)
        self.include_type_labels_label = tk.Label(self.options_frame, text='include_type_labels', font=font_settings)
        self.include_type_labels_label.grid(column=0, row=2)
        
        # neighbor_charge_constraint drop down menu
        styles = ['check-neighbors', 'accumulate-carbon', 'accumulate-pi-electron', 'accumulate-neighbor', 'none']
        self.neighbor_charge_constraint = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/2.75), font=font_settings)
        self.neighbor_charge_constraint.current(styles.index(neighbor_charge_constraint))
        self.neighbor_charge_constraint.grid(column=1, row=3, columnspan=2)
        self.neighbor_charge_constraint_label = tk.Label(self.options_frame, text='neighbor_charge_constraint', font=font_settings)
        self.neighbor_charge_constraint_label.grid(column=1, row=2, columnspan=2)
        
        # types2convert entry
        self.types2convert = tk.Entry(self.options_frame, width=int(maxwidth/1.125), font=font_settings)
        self.types2convert.insert(0, ','.join([str(i) for i in types2convert]))
        self.types2convert.grid(column=3, row=3, columnspan=4)
        self.types2convert_label = tk.Label(self.options_frame, text='types2convert (comma separated w/no whitespace)', font=font_settings)
        self.types2convert_label.grid(column=3, row=2, columnspan=4)
        

        
        # Add padding to all frames in self.inputs_frame
        for widget in self.options_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))

        #------------#
        # Run button #
        #------------#
        self.run = tk.Button(self.frame, text='Run LUNAR/add_pi_electrons.py', font=font_settings, command=self.run_LUNAR)
        self.run.grid(row=4, column=0, columnspan=2, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
        #-----------------#
        # update defaults #
        #-----------------#
        self.update = tk.Button(self.frame, text='Save the current GUI settings as the default GUI settings', font=font_settings, command=self.update_py_script)
        self.update.grid(row=5, column=0, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
        #------------#
        # Quick help #
        #------------#
        self.quick_help = tk.Button(self.frame, text='Quick help', font=font_settings, command=self.quickhelp)
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
            txt = os.path.join(self.pwd, 'src/GUI_quick_help_pages/add_pi_electrons.txt')
            logged = []
            with open(txt, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    # Strip comment's and split by whitespace
                    line = line.split('#')[0]
                    line = line.rstrip()
                    logged.append(line)
        except: # except something failed
            logged.append('FAILED to read LUNAR/src/GUI_quick_help_pages/add_pi_electrons.txt document.')
            logged.append('Most likely cause is the add_pi_electrons.txt file was renamed in the')
            logged.append('LUNAR/src/GUI_quick_help_pages/ directory or directory names were changed.')
        self.popup(logged, title='Quick help')
        return  
    
    # Function to get filepath for topofile
    def topofile_path(self):
        ftypes = (('data files', '*.data'), ('all files', '*.*'))
        path = filedialog.askopenfilename(title='Open topofile?', filetypes=ftypes)
        if path:
            path = os.path.relpath(path)
            self.topofile.delete(0, tk.END); self.topofile.insert(0, path);
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
        print('Terminating add_pi_electrons GUI'); self.root.destroy();
        return
    
    # Function to run imported LUNAR functions
    def run_LUNAR(self): 
        # Set up log
        log = io_functions.LUNAR_logger()
        valid_inputs = True
        
        # Get information from GUI
        boolean = {'False':False, 'True':True}
        topofile = self.topofile.get()
        parent_directory = self.parent_directory.get() 
        newfile = self.newfile.get()
        atom_style = self.atom_style.get()
        try: 
            types2convert = s2s_i2i(self.types2convert.get().split(','))
            if types2convert == ['']:
                log.GUI_error('ERROR types2convert is an empty list')
                valid_inputs = False
        except: 
            types2convert = []
            log.GUI_error('ERROR types2convert is an empty list')
            valid_inputs = False
        convert2cg1 = boolean[self.convert2cg1.get()]
        reset_charges = boolean[self.reset_charges.get()]
        net_zero_charge = boolean[self.net_zero_charge.get()]
        add_pi_electrons = boolean[self.add_pi_electrons.get()]
        reset_simulation_cell = boolean[self.reset_simulation_cell.get()]
        include_type_labels = boolean[self.include_type_labels.get()]
        neighbor_charge_constraint = self.neighbor_charge_constraint.get()
        
        # Run LUNAR/add_pi_electrons
        if valid_inputs:
            try: 
                inputs = (topofile, types2convert, atom_style, reset_charges, net_zero_charge, convert2cg1, add_pi_electrons,
                          parent_directory, newfile, include_type_labels, neighbor_charge_constraint, reset_simulation_cell, [], log)
                
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
        parent_directory = io_functions.path_to_string(self.parent_directory.get()) 
        newfile = self.newfile.get()
        atom_style = self.atom_style.get()
        try: types2convert = s2s_i2i(self.types2convert.get().split(','))
        except: types2convert = []
        convert2cg1 = boolean[self.convert2cg1.get()]
        reset_charges = boolean[self.reset_charges.get()]
        net_zero_charge = boolean[self.net_zero_charge.get()]
        add_pi_electrons = boolean[self.add_pi_electrons.get()]
        reset_simulation_cell = boolean[self.reset_simulation_cell.get()]
        include_type_labels = boolean[self.include_type_labels.get()]
        neighbor_charge_constraint = self.neighbor_charge_constraint.get()
        
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
                if line.startswith('parent_directory') and inputsflag:
                    line = psm.parse_and_modify(line, parent_directory, stringflag=True, splitchar='=')
                if line.startswith('include_type_labels') and inputsflag:
                    line = psm.parse_and_modify(line, include_type_labels, stringflag=False, splitchar='=')
                if line.startswith('reset_simulation_cell') and inputsflag:
                    line = psm.parse_and_modify(line, reset_simulation_cell, stringflag=False, splitchar='=')
                if line.startswith('newfile') and inputsflag:
                    line = psm.parse_and_modify(line, newfile, stringflag=True, splitchar='=')
                if line.startswith('atom_style') and inputsflag:
                    line = psm.parse_and_modify(line, atom_style, stringflag=True, splitchar='=')
                if line.startswith('neighbor_charge_constraint') and inputsflag:
                    line = psm.parse_and_modify(line, neighbor_charge_constraint, stringflag=True, splitchar='=')
                if line.startswith('types2convert') and inputsflag:
                    line = psm.parse_and_modify(line, str(types2convert), stringflag=False, splitchar='=')
                if line.startswith('convert2cg1') and inputsflag:
                    line = psm.parse_and_modify(line, convert2cg1, stringflag=False, splitchar='=')
                if line.startswith('reset_charges') and inputsflag:
                    line = psm.parse_and_modify(line, reset_charges, stringflag=False, splitchar='=')
                if line.startswith('net_zero_charge') and inputsflag:
                    line = psm.parse_and_modify(line, net_zero_charge, stringflag=False, splitchar='=')
                if line.startswith('add_pi_electrons') and inputsflag:
                    line = psm.parse_and_modify(line, add_pi_electrons, stringflag=False, splitchar='=')
                if line.startswith('if __name__ == "__main__":'): inputsflag = False
                f.write(line)
        return