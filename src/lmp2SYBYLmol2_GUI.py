# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
May 16, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.GUI_scale_settings as GUI_scale_settings
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
import sys
import os


# setting path to files called by other codes in src
sys.path.append('../')
from lmp2SYBYLmol2 import main


###########################
# LUNAR/lmp2SYBYLmol2 GUI #
###########################
class GUI:
    def __init__(self, topofile, parent_directory, remove_PBC_bonds, mass_map, addbox, GUI_zoom):
        self.mass_map = mass_map
        
        # Find present working directory
        self.pwd = os.getcwd()
        self.filename = os.path.join(self.pwd, 'lmp2SYBYLmol2.py')

        # Initialize window
        self.root = tk.Tk()
        self.root.title('LUNAR/lmp2SYBYLmol2.py GUI v1.0')
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
        self.topofile = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.topofile.insert(0, topofile)
        self.topofile.grid(column=1, row=0, columnspan=3)
        self.topofile_button = tk.Button(self.inputs_frame, text='topofile', font=font_settings, command=self.topofile_path)
        self.topofile_button.grid(column=0, row=0)
        
        # parent_directory entry
        self.parent_directory = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.parent_directory.insert(0, parent_directory)
        self.parent_directory.grid(column=1, row=1, columnspan=3)
        self.dir_button = tk.Button(self.inputs_frame, text='parent_directory', font=font_settings, command=self.directory_path)
        self.dir_button.grid(column=0, row=1)
        
        # remove_PBC_bonds drop down menu
        styles = [True, False]
        self.remove_PBC_bonds = ttk.Combobox(self.inputs_frame, values=styles, width=int(maxwidth/3), font=font_settings)
        self.remove_PBC_bonds.current(styles.index(remove_PBC_bonds))
        self.remove_PBC_bonds.grid(column=1, row=2)
        self.remove_PBC_bonds_label = tk.Label(self.inputs_frame, text='remove_PBC_bonds', font=font_settings)
        self.remove_PBC_bonds_label.grid(column=0, row=2)
        
        # remove_PBC_bonds drop down menu
        styles = [True, False]
        self.addbox = ttk.Combobox(self.inputs_frame, values=styles, width=int(maxwidth/3), font=font_settings)
        self.addbox.current(styles.index(addbox))
        self.addbox.grid(column=3, row=2)
        self.addbox_label = tk.Label(self.inputs_frame, text='addbox', font=font_settings)
        self.addbox_label.grid(column=2, row=2)
        

        # Add padding to all frames in self.inputs_frame
        for widget in self.inputs_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=ypadding)

        #------------#
        # Run button #
        #------------#
        self.run = tk.Button(self.frame, text='Run LUNAR/lmp2SYBYLmol2.py', font=font_settings, command=self.run_LUNAR)
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
            txt = os.path.join(self.pwd, 'src/GUI_quick_help_pages/lmp2SYBYLmol2.txt')
            logged = []
            with open(txt, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    # Strip comment's and split by whitespace
                    line = line.split('#')[0]
                    line = line.rstrip()
                    logged.append(line)
        except: # except something failed
            logged.append('FAILED to read LUNAR/src/GUI_quick_help_pages/lmp2SYBYLmol2.txt document.')
            logged.append('Most likely cause is the lmp2SYBYLmol2.txt file was renamed in the')
            logged.append('LUNAR/src/GUI_quick_help_pages/ directory or directory names were changed.')
        self.popup(logged, title='Quick help')
        return   
    
    # Function to get filepath for topofile
    def topofile_path(self):
        ftypes = (('data files', '*.data *.data.gz'), ('all files', '*.*'))
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
        print('Terminating lmp2SYBYLmol2 GUI'); self.root.destroy();
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
        remove_PBC_bonds = boolean[self.remove_PBC_bonds.get()]
        addbox = boolean[self.addbox.get()]
        mass_map = self.mass_map
        
        # Run LUNAR/lmp2SYBYLmol2
        if valid_inputs:
            try: 
                inputs = (topofile, parent_directory, remove_PBC_bonds, mass_map, addbox, log)
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
        remove_PBC_bonds = boolean[self.remove_PBC_bonds.get()]
        addbox = boolean[self.addbox.get()]
        
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
                if line.startswith('remove_PBC_bonds') and inputsflag:
                    line = psm.parse_and_modify(line, remove_PBC_bonds, stringflag=False, splitchar='=')
                if line.startswith('addbox') and inputsflag:
                    line = psm.parse_and_modify(line, addbox, stringflag=False, splitchar='=')
                if line.startswith('if __name__ == "__main__":') or line.startswith('def main('): inputsflag = False
                f.write(line)
        return