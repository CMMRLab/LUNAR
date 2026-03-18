# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
October 21, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.GUI_scale_settings as GUI_scale_settings
import src.io_functions as io_functions
import src.array_lmp_script.mode_from_script as mode_from_script
import src.array_lmp_script.main as main
from tkinter.scrolledtext import ScrolledText
from tkinter import filedialog
from tkinter import Toplevel
from tkinter import ttk
import tkinter as tk
import traceback
import math
import os



###############
# plotter GUI #
###############
class GUI:
    def __init__(self, settings, GUI_zoom):
        
        # Configure log (default is level='production', switch to 'debug' if debuging)
        self.log=io_functions.LUNAR_logger()
        self.log.configure(level='production')

        # Find present working directory
        self.pwd = os.getcwd()
        self.modespath = settings['modes-dir']
        
        # Set defaults
        self.settings = settings
        
        try:
            module = main.import_file(settings['mode'])
            self.mode = module.mode
        except:            
            self.mode = {'lmp_script': 'UPDATE-ME',
                         'batch_script': '$<lmp_script>.sh',
                         'auto_script': 'auto_submit_$<lmp_script>',
                         'build_script': '$<lmp_script>',
                         'build_directory': 'UPDATE-ME',
                         'lmp_script_find': '$<lmp_script>',
                         'build_script_find': '$<build_script>',
                         'auto_submit': 'qsub',
                         'lmp_array': {}, 
                         'batch_array': {},
                         'copy_files': []}

            self.log.GUI_error(f'ERROR loading mode file {settings["mode"]}. Internally deriving default settings. Likely launching from outside of LUNAR directory.')

        
        # Initialize window
        self.root = tk.Tk()
        self.root.title('LUNAR/array_lmp_script.py GUI v1.0')
        self.root.resizable(width=False, height=False)
        #self.root.geometry('600x400')
        
        # Setup tabs
        self.tabControl = ttk.Notebook(self.root)
        self.tab1 = ttk.Frame(self.tabControl)
        self.tab2 = ttk.Frame(self.tabControl)
        self.tabControl.add(self.tab1, text='Main')
        self.tabControl.add(self.tab2, text='Additional Settings')
        self.tabControl.pack(expand=1, fill='both')
        #self.root.geometry('600x400')
        
        # Initalize main frame
        self.frame1 = tk.Frame(self.tab1)
        self.frame1.pack()
        self.frame2 = tk.Frame(self.tab2)
        self.frame2.pack()
        
        # Initalize main frame
        #self.frame = tk.Frame(self.root)
        #self.frame.pack()
        
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
        self.maxwidth = 110

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
        self.xpadding = int(math.ceil(GUI_SF*self.xpadding))
        self.ypadding = int(math.ceil(GUI_SF*self.ypadding))
        self.maxwidth = int(math.ceil(GUI_SF*self.maxwidth))
        self.font_settings = (self.font_type, font_size)
        self.GUI_zoom = GUI_zoom
        
        
        #--------------#
        # Inputs frame #
        #--------------#
        # Initalize  inputs frame
        self.inputs_frame = tk.LabelFrame(self.frame1, text='Inputs', font=self.font_settings)
        self.inputs_frame.grid(row=0, column=0, padx=self.xpadding, pady=self.ypadding)
        
        # lmp_script selection button
        self.lmp_script = tk.Entry(self.inputs_frame, width=int(1*self.maxwidth), font=self.font_settings)
        self.lmp_script.insert(0, self.mode['lmp_script'])
        self.lmp_script.grid(column=1, row=0, columnspan=1)
        self.lmp_script_button = tk.Button(self.inputs_frame, text='lmp_script', font=self.font_settings, command=self.lmp_script_path)
        self.lmp_script_button.grid(column=0, row=0)
        
        # lmp_script_find
        self.lmp_script_find= tk.Entry(self.inputs_frame, width=int(0.375*self.maxwidth), font=self.font_settings)
        self.lmp_script_find.insert(0, self.mode['lmp_script_find'])
        self.lmp_script_find.grid(column=4, row=0)
        self.lmp_script_find_label = tk.Label(self.inputs_frame, text='lmp_script_find (e.g. $<lmp_script>)', font=self.font_settings)
        self.lmp_script_find_label.grid(column=3, row=0)
        
        # batch_script selection button
        self.batch_script = tk.Entry(self.inputs_frame, width=int(1*self.maxwidth), font=self.font_settings)
        self.batch_script.insert(0, self.mode['batch_script'])
        self.batch_script.grid(column=1, row=1, columnspan=1)
        self.batch_script_button = tk.Button(self.inputs_frame, text='batch_script', font=self.font_settings, command=self.batch_script_path)
        self.batch_script_button.grid(column=0, row=1)
        
        # modes
        self.modefile = tk.Entry(self.inputs_frame, width=int(1*self.maxwidth), font=self.font_settings)
        self.modefile.insert(0, settings['mode'])
        self.modefile.grid(column=1, row=2)
        self.modefile_button = tk.Button(self.inputs_frame, text='mode file', font=self.font_settings, command=self.modefile_path)
        self.modefile_button.grid(column=0, row=2)        
        
        # load_replace_logfile drop down menu
        self.bools_list = [True, False]
        self.load_replace_logfile = ttk.Combobox(self.inputs_frame, values=self.bools_list, width=int(0.3525*self.maxwidth), font=self.font_settings)
        self.load_replace_logfile.current(self.bools_list.index(settings['replace_lmp_script_when_loading_mode']))
        self.load_replace_logfile.grid(column=4, row=2)
        self.load_replace_logfile_label = tk.Label(self.inputs_frame, text='Replace lmp_script when loading mode', font=self.font_settings)
        self.load_replace_logfile_label.grid(column=3, row=2)
        
        # auto_script
        self.auto_script = tk.Entry(self.inputs_frame, width=int(1*self.maxwidth), font=self.font_settings)
        self.auto_script.insert(0, self.mode['auto_script'])
        self.auto_script.grid(column=1, row=3, columnspan=1)
        self.auto_script_label = tk.Label(self.inputs_frame, text='auto_script', font=self.font_settings)
        self.auto_script_label.grid(column=0, row=3)
        
        # auto_submit
        self.auto_submit = tk.Entry(self.inputs_frame, width=int(0.375*self.maxwidth), font=self.font_settings)
        self.auto_submit.insert(0, self.mode['auto_submit'])
        self.auto_submit.grid(column=4, row=3)
        self.auto_submit_label = tk.Label(self.inputs_frame, text='auto_submit (e.g. qsub)', font=self.font_settings)
        self.auto_submit_label.grid(column=3, row=3)
        
        # build_script
        self.build_script = tk.Entry(self.inputs_frame, width=int(1*self.maxwidth), font=self.font_settings)
        self.build_script.insert(0, self.mode['build_script'])
        self.build_script.grid(column=1, row=4, columnspan=1)
        self.build_script_label = tk.Label(self.inputs_frame, text='build_script', font=self.font_settings)
        self.build_script_label.grid(column=0, row=4)
        
        # build_script_find
        self.build_script_find= tk.Entry(self.inputs_frame, width=int(0.375*self.maxwidth), font=self.font_settings)
        self.build_script_find.insert(0, self.mode['build_script_find'])
        self.build_script_find.grid(column=4, row=4)
        self.build_script_find_label = tk.Label(self.inputs_frame, text='build_script_find (e.g. $<build_script>)', font=self.font_settings)
        self.build_script_find_label.grid(column=3, row=4)
        
        # build_directory
        self.build_directory = tk.Entry(self.inputs_frame, width=int(1.75*self.maxwidth), font=self.font_settings)
        self.build_directory.insert(0, self.mode['build_directory'])
        self.build_directory.grid(column=1, row=5, columnspan=4)
        self.build_directory_label = tk.Label(self.inputs_frame, text='build_directory', font=self.font_settings)
        self.build_directory_label.grid(column=0, row=5)

        # Add padding to all frames in self.inputs_frame
        for widget in self.inputs_frame.winfo_children():
            widget.grid_configure(padx=int(self.xpadding/2), pady=int(self.ypadding/3))
            
            
        #-----------------#
        # lmp_array frame #
        #-----------------#
        # Initalize lmp_array frame
        self.lmp_array_frame = tk.LabelFrame(self.frame1, text='lmp_array', font=self.font_settings)
        self.lmp_array_frame.grid(row=1, column=0, sticky='news', padx=self.xpadding, pady=self.ypadding)
        
        # file selection labels
        self.lmp_array_find_label = tk.Label(self.lmp_array_frame, text='Find (e.g. $<replicate> or $<seed> or $<ints> or ...)', font=self.font_settings)
        self.lmp_array_find_label.grid(column=0, row=0)
        
        self.lmp_array_replace_label = tk.Label(self.lmp_array_frame, text='Replace (e.g. [1, 2, 3] or [12345] or range(1, 5) or ...)', font=self.font_settings)
        self.lmp_array_replace_label.grid(column=1, row=0, columnspan=2)
        
        
        # file selection button and qty
        self.nlmp_array = len(self.mode['lmp_array'])
        self.lmp_array_find = []; self.lmp_array_replace = []
        for n, (key, value) in enumerate(self.mode['lmp_array'].items()):

            self.lmp_find = tk.Entry(self.lmp_array_frame, width=int(0.65*self.maxwidth), font=self.font_settings)
            self.lmp_find.grid(column=0, row=n+1)
            self.lmp_find.insert(0, key)
            self.lmp_array_find.append(self.lmp_find)
            
            self.lmp_replace = tk.Entry(self.lmp_array_frame, width=int(1.25*self.maxwidth), font=self.font_settings)
            self.lmp_replace.grid(column=1, row=n+1, columnspan=2)
            self.lmp_replace.insert(0, str(value))
            self.lmp_array_replace.append(self.lmp_replace)
            
        # Button to add a file
        self.add_button_lmp = tk.Button(self.lmp_array_frame, text='add find/replace to stack', font=self.font_settings, width=int(self.maxwidth/1.75), command=self.add2stack_lmp)
        self.add_button_lmp.grid(column=0, row=self.nlmp_array+1, sticky='news', columnspan=1)
            
        # Button to remove a file
        self.remove_button_lmp = tk.Button(self.lmp_array_frame, text='remove last find/replace from stack', font=self.font_settings, width=int(self.maxwidth/3.5), command=self.remove_last_lmp)
        self.remove_button_lmp.grid(column=1, row=self.nlmp_array+1, sticky='news', columnspan=1)
        
        # Button to clear all files
        self.clear_button_lmp = tk.Button(self.lmp_array_frame, text='clear stack', font=self.font_settings, width=int(self.maxwidth/3.5), command=self.clear_all_lmp)
        self.clear_button_lmp.grid(column=2, row=self.nlmp_array+1, sticky='news', columnspan=1)
        
        # Add padding to all frames in self.lmp_array_frame
        for widget in self.lmp_array_frame.winfo_children():
            widget.grid_configure(padx=int(self.xpadding/4), pady=int(self.ypadding/3))
            
            
        #-------------------#
        # batch_array frame #
        #-------------------#
        # Initalize batch_array frame
        self.batch_array_frame = tk.LabelFrame(self.frame2, text='batch_array', font=self.font_settings)
        self.batch_array_frame.grid(row=2, column=0, sticky='news', padx=self.xpadding, pady=self.ypadding)
        
        # file selection labels
        self.batch_array_find_label = tk.Label(self.batch_array_frame, text='Find (e.g. $<ncores> or $<que> or $<input> or $<output> or ...)', font=self.font_settings)
        self.batch_array_find_label.grid(column=0, row=0)
        
        self.batch_array_replace_label = tk.Label(self.batch_array_frame, text='Replace (e.g. 16 or long.q or $<lmp_scipt> $<lmp_scipt>.out or ...)', font=self.font_settings)
        self.batch_array_replace_label.grid(column=1, row=0, columnspan=2)
        
        
        # file selection button and qty
        self.nbatch_array = len(self.mode['batch_array'])
        self.batch_array_find = []; self.batch_array_replace = []
        for n, (key, value) in enumerate(self.mode['batch_array'].items()):

            self.batch_find = tk.Entry(self.batch_array_frame, width=int(0.65*self.maxwidth), font=self.font_settings)
            self.batch_find.grid(column=0, row=n+1)
            self.batch_find.insert(0, key)
            self.batch_array_find.append(self.batch_find)
            
            self.batch_replace = tk.Entry(self.batch_array_frame, width=int(1.25*self.maxwidth), font=self.font_settings)
            self.batch_replace.grid(column=1, row=n+1, columnspan=2)
            self.batch_replace.insert(0, str(value))
            self.batch_array_replace.append(self.batch_replace)
            
        # Button to add a file
        self.add_button_batch = tk.Button(self.batch_array_frame, text='add find/replace to stack', font=self.font_settings, width=int(self.maxwidth/1.75), command=self.add2stack_batch)
        self.add_button_batch.grid(column=0, row=self.nbatch_array+1, sticky='news', columnspan=1)
            
        # Button to remove a file
        self.remove_button_batch = tk.Button(self.batch_array_frame, text='remove last find/replace from stack', font=self.font_settings, width=int(self.maxwidth/3.5), command=self.remove_last_batch)
        self.remove_button_batch.grid(column=1, row=self.nbatch_array+1, sticky='news', columnspan=1)
        
        # Button to clear all files
        self.clear_button_batch = tk.Button(self.batch_array_frame, text='clear stack', font=self.font_settings, width=int(self.maxwidth/3.5), command=self.clear_all_batch)
        self.clear_button_batch.grid(column=2, row=self.nbatch_array+1, sticky='news', columnspan=1)
        
        # Add padding to all frames in self.batch_array_frame
        for widget in self.batch_array_frame.winfo_children():
            widget.grid_configure(padx=int(self.xpadding/4), pady=int(self.ypadding/3))
            
        #------------------#
        # copy_files frame #
        #------------------#
        # Initalize copy_files frame
        self.copy_files_frame = tk.LabelFrame(self.frame2, text='copy_files with find/replace operations (e.g. MyData.data or MyMolecule_replicate_$<replicate> or *.data or * ...)', font=self.font_settings)
        self.copy_files_frame.grid(row=3, column=0, sticky='news', padx=self.xpadding, pady=self.ypadding)
        
        # file selection button and qty
        self.ncopy_files = len(self.mode['copy_files'])
        self.copy_files = []
        for n, file in enumerate(self.mode['copy_files']):
            self.copy_file = tk.Entry(self.copy_files_frame, width=int(1.9*self.maxwidth), font=self.font_settings)
            self.copy_file.grid(column=0, row=n+1, columnspan=3)
            self.copy_file.insert(0, file)
            self.copy_files.append(self.copy_file)
            
        # Button to add a file
        self.add_button_copy = tk.Button(self.copy_files_frame, text='add file(s) to stack', font=self.font_settings, width=int(self.maxwidth/3), command=self.add2stack_copy)
        self.add_button_copy.grid(column=0, row=self.ncopy_files+1, sticky='news', columnspan=1)
            
        # Button to remove a file
        self.remove_button_copy = tk.Button(self.copy_files_frame, text='remove last file stack', font=self.font_settings, width=int(self.maxwidth/3), command=self.remove_last_copy)
        self.remove_button_copy.grid(column=1, row=self.ncopy_files+1, sticky='news', columnspan=1)
        
        # Button to clear all files
        self.clear_button_copy = tk.Button(self.copy_files_frame, text='clear stack', font=self.font_settings, width=int(self.maxwidth/3), command=self.clear_all_copy)
        self.clear_button_copy.grid(column=2, row=self.ncopy_files+1, sticky='news', columnspan=1)
        
        # Add padding to all frames in self.copy_files_frame
        for widget in self.copy_files_frame.winfo_children():
            widget.grid_configure(padx=int(self.xpadding/4), pady=int(self.ypadding/3))

            
        #---------------------------#
        # Initalize execute frame:1 #
        #---------------------------#
        self.execute_frame1 = tk.LabelFrame(self.frame1, borderwidth=0, highlightthickness=0)
        self.execute_frame1.grid(column=0, row=4, padx=self.xpadding, pady=self.ypadding)
        
        # run button
        self.run_btn = tk.Button(self.execute_frame1, width=int(1.9*self.maxwidth), text='Run LUNAR/array_lmp_script.py', command=self.run_main, font=self.font_settings)
        self.run_btn.grid(row=0, column=0, columnspan=2, sticky='news')
        
        # save_mode
        self.save_mode_btn = tk.Button(self.execute_frame1, width=int(self.maxwidth/1.75), text='save settings as mode', command=self.save_mode, font=self.font_settings)
        self.save_mode_btn.grid(row=1, column=0, sticky='news')
        
        # Quick help button
        self.quick_help = tk.Button(self.execute_frame1, width=int(self.maxwidth/1.75), text='Quick help', font=self.font_settings, command=self.quickhelp)
        self.quick_help.grid(row=1, column=1, sticky='news')

        
        # Add padding to all frames in self.execute_frame
        for widget in self.execute_frame1.winfo_children():
            widget.grid_configure(padx=int(self.xpadding/2), pady=int(self.ypadding/3))
            
        #---------------------------#
        # Initalize execute frame:2 #
        #---------------------------#
        self.execute_frame2 = tk.LabelFrame(self.frame2, borderwidth=0, highlightthickness=0)
        self.execute_frame2.grid(column=0, row=4, padx=self.xpadding, pady=self.ypadding)
        
        # run button
        self.run_btn = tk.Button(self.execute_frame2, width=int(1.9*self.maxwidth), text='Run LUNAR/array_lmp_script.py', command=self.run_main, font=self.font_settings)
        self.run_btn.grid(row=0, column=0, columnspan=2, sticky='news')
        
        # save_mode
        self.save_mode_btn = tk.Button(self.execute_frame2, width=int(self.maxwidth/1.75), text='save settings as mode', command=self.save_mode, font=self.font_settings)
        self.save_mode_btn.grid(row=1, column=0, sticky='news')
        
        # Quick help button
        self.quick_help = tk.Button(self.execute_frame2, width=int(self.maxwidth/1.75), text='Quick help', font=self.font_settings, command=self.quickhelp)
        self.quick_help.grid(row=1, column=1, sticky='news')

        
        # Add padding to all frames in self.execute_frame
        for widget in self.execute_frame2.winfo_children():
            widget.grid_configure(padx=int(self.xpadding/2), pady=int(self.ypadding/3))
        
        
        #------------------------#
        # Run mainloop and close #
        #------------------------#
        self.root.protocol('WM_DELETE_WINDOW', self.closing)
        self.root.mainloop()
        
    
    ################################
    # Functions to call as methods #
    ################################             
    # Function to construct mode from GUI settings
    def construct_mode_dict(self):
        mode = {}
        mode['lmp_script'] = io_functions.path_to_string(self.lmp_script.get())
        mode['batch_script'] = io_functions.path_to_string(self.batch_script.get())
        mode['auto_script'] = io_functions.path_to_string(self.auto_script.get())
        mode['build_script'] = io_functions.path_to_string(self.build_script.get())
        mode['build_directory'] = io_functions.path_to_string(self.build_directory.get())
        mode['auto_submit'] = io_functions.path_to_string(self.auto_submit.get())
        mode['lmp_script_find'] = self.lmp_script_find.get()
        mode['build_script_find'] = self.build_script_find.get()
        mode['lmp_array'] = {}
        mode['batch_array'] = {}
        mode['copy_files'] = []
        for i, j  in zip(self.lmp_array_find, self.lmp_array_replace):
            find = i.get(); replace = j.get()
            if find != '' and replace != '': 
                mode['lmp_array'][find] = replace
                
        for i, j  in zip(self.batch_array_find, self.batch_array_replace):
            find = i.get(); replace = j.get()
            if find != '' and replace != '': 
                mode['batch_array'][find] = replace
                
        for i in self.copy_files:
            file = i.get()
            if file != '':
                mode['copy_files'].append(file)
                
        # for n, i in enumerate(mode, 1):
        #     j = mode[i]
        #     print('{} "{}" = "{}"'.format(n, i, str(j)))
        return mode
    
    
    # Quick help button
    def quickhelp(self):
        try: # Try to get text from GUI_help_page.txt file
            txt = os.path.join(self.pwd, 'src/GUI_quick_help_pages/array_lmp_script.txt')
            logged = []
            with open(txt, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    # Strip comment's if '!#' not in line,
                    # else use the '!#' combination to keep
                    # the '#' character
                    if '!#' not in line:
                        line = line.split('#')[0]
                    else:
                        line = line.replace('!#', '#')
                    line = line.rstrip()
                    logged.append(line)
        except: # except something failed
            logged.append('FAILED to read LUNAR/src/GUI_quick_help_pages/array_lmp_script.txt document.')
            logged.append('Most likely cause is the array_lmp_script.txt file was renamed in the')
            logged.append('LUNAR/src/GUI_quick_help_pages/ directory or directory names were changed.')
        self.popup(logged, title='Quick help')
        return  
    
    # Function to add to analysis options
    def add2stack_copy(self):
        ftypes = (('all files', '*.*'), ('data files', '*.data *.data.gz'))
        paths = filedialog.askopenfilenames(title='Open file(s)?', filetypes=ftypes)
        if paths:
            for path in paths:
                file = os.path.relpath(path)
                file = io_functions.path_to_string(file)
                try:
                    lmp_script = io_functions.path_to_string(self.lmp_script.get())
                    lmp_directory = os.path.dirname(lmp_script)
                    
                    relative = os.path.relpath(file, start=lmp_directory)
                    relative = io_functions.path_to_string(relative)
                except: relative = file
                
                # Any open slots
                any_available = False
                n_available = 0
                for n, copy_file in enumerate(self.copy_files): 
                    if copy_file.get() == '':
                        any_available = True
                        n_available = n
                        break
                if any_available:
                    self.copy_files[n_available].insert(0, relative)
                else:
                    try: self.add_overloaded_analysis_copy(file=relative)
                    except: print(f'GUI failed to add a copy file to the stack: {file}')
        return
    
    # Function to add files to GUI, during overload conditions
    def add_overloaded_analysis_copy(self, file=''):    
        # adjust based on GUI_SF
        GUI_SF = self.GUI_zoom/100
        xpadding = int(math.ceil(GUI_SF*self.xpadding))
        ypadding = int(math.ceil(GUI_SF*self.ypadding))
        
        # Add file box
        self.ncopy_files += 1
        
        self.copy_file = tk.Entry(self.copy_files_frame, width=int(1.9*self.maxwidth), font=self.font_settings)
        self.copy_file.grid(column=0, row=self.ncopy_files, columnspan=3)
        self.copy_file.insert(0, file)
        self.copy_files.append(self.copy_file)
        
        # adjust packing of other things in inputs frame
        self.add_button_copy.grid(column=0, row=self.ncopy_files+1, sticky='news', columnspan=1)
        self.remove_button_copy.grid(column=1, row=self.ncopy_files+1, sticky='news', columnspan=1)
        self.clear_button_copy.grid(column=2, row=self.ncopy_files+1, sticky='news', columnspan=1)
        
        # Add padding to all frames in self.copy_files_frame
        for widget in self.copy_files_frame.winfo_children():
            widget.grid_configure(padx=int(xpadding/4), pady=int(ypadding/3))
        return
    
    # Function to remove last anaylsis
    def remove_last_copy(self):
        used = []
        for index, i in enumerate(self.copy_files):
            file = i.get()
            if file != '': used.append(index)
        if used:
            last_add = max(used)
            self.copy_files[last_add].delete(0, tk.END)
        else: print('No find or replace left to remove')
        return
    
    # Function to clear all files
    def clear_all_copy(self):
        for n, i in enumerate(self.copy_files):
            try: self.copy_files[n].delete(0, tk.END)
            except: pass
        
        self.destroy_copy_rows(keep_first_row=True)
        return
    
    # Method to destroy copy files rows
    def destroy_copy_rows(self, keep_first_row=False):
        """
        Destroys the dynamic peak row widgets (Entries/Comboboxes) and clears tracking lists.
        If keep_first_row=True, keeps row 1 (index 0) widgets and destroys the rest.
        """
        if keep_first_row: start = 1
        else: start = 0
    
        # destroy widgets from the end (safer)
        for idx in range(len(self.copy_files)-1, start-1, -1):
            w = self.copy_files[idx]
            try: w.destroy()
            except: pass
    
        # shrink lists and update n-to-match
        self.copy_files  = self.copy_files[:start]
        self.ncopy_files = len(self.copy_files)
        return
    
    # Function to add to analysis options
    def add2stack_batch(self):
        try: self.add_overloaded_analysis_batch()
        except: print('GUI failed to add additionaly analysis to stack')
        return
    
    # Function to add files to GUI, during overload conditions
    def add_overloaded_analysis_batch(self):    
        # adjust based on GUI_SF
        GUI_SF = self.GUI_zoom/100
        xpadding = int(math.ceil(GUI_SF*self.xpadding))
        ypadding = int(math.ceil(GUI_SF*self.ypadding))
        
        # Add file box
        self.nbatch_array += 1
        
        self.batch_find = tk.Entry(self.batch_array_frame, width=int(0.65*self.maxwidth), font=self.font_settings)
        self.batch_find.grid(column=0, row=self.nbatch_array)
        self.batch_find.insert(0, '')
        self.batch_array_find.append(self.batch_find)
        
        self.batch_replace = tk.Entry(self.batch_array_frame, width=int(1.25*self.maxwidth), font=self.font_settings)
        self.batch_replace.grid(column=1, row=self.nbatch_array, columnspan=2)
        self.batch_replace.insert(0, '')
        self.batch_array_replace.append(self.batch_replace)
        
        # adjust packing of other things in inputs frame
        self.add_button_batch.grid(column=0, row=self.nbatch_array+1, sticky='news', columnspan=1)
        self.remove_button_batch.grid(column=1, row=self.nbatch_array+1, sticky='news', columnspan=1)
        self.clear_button_batch.grid(column=2, row=self.nbatch_array+1, sticky='news', columnspan=1)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.batch_array_frame.winfo_children():
            widget.grid_configure(padx=int(xpadding/4), pady=int(ypadding/3))
        return
    
    # Function to remove last anaylsis
    def remove_last_batch(self):
        used = []
        for index, (i, j)  in enumerate(zip(self.batch_array_find, self.batch_array_replace)):
            find = i.get(); replace = j.get()
            if find != '' or replace != '': used.append(index)
        if used:
            last_add = max(used)
            self.batch_array_find[last_add].delete(0, tk.END)
            self.batch_array_replace[last_add].delete(0, tk.END)
        else: print('No find or replace left to remove')
        return
    
    # Function to clear all files
    def clear_all_batch(self):
        for n, i in enumerate(self.batch_array_find):
            try: self.batch_array_find[n].delete(0, tk.END)
            except: pass
        for n, i in enumerate(self.batch_array_replace):
            try: self.batch_array_replace[n].delete(0, tk.END)
            except: pass  
        
        self.destroy_batch_rows(keep_first_row=True)
        return
    
    # Method to destroy batch rows
    def destroy_batch_rows(self, keep_first_row=False):
        """
        Destroys the dynamic peak row widgets (Entries/Comboboxes) and clears tracking lists.
        If keep_first_row=True, keeps row 1 (index 0) widgets and destroys the rest.
        """
        if keep_first_row: start = 1
        else: start = 0
    
        # destroy widgets from the end (safer)
        for idx in range(len(self.batch_array_find)-1, start-1, -1):
            for w in (self.batch_array_find[idx], self.batch_array_replace[idx]):
                try: w.destroy()
                except: pass
    
        # shrink lists and update n-to-match
        self.batch_array_find    = self.batch_array_find[:start]
        self.batch_array_replace = self.batch_array_replace[:start]
        self.nbatch_array        = len(self.batch_array_find)
        return
    
    # Function to add to analysis options
    def add2stack_lmp(self):
        try: self.add_overloaded_analysis_lmp()
        except: print('GUI failed to add additionaly analysis to stack')
        return
    
    # Function to add files to GUI, during overload conditions
    def add_overloaded_analysis_lmp(self):    
        # adjust based on GUI_SF
        GUI_SF = self.GUI_zoom/100
        xpadding = int(math.ceil(GUI_SF*self.xpadding))
        ypadding = int(math.ceil(GUI_SF*self.ypadding))
        
        # Add file box
        self.nlmp_array += 1
        
        self.lmp_find = tk.Entry(self.lmp_array_frame, width=int(0.65*self.maxwidth), font=self.font_settings)
        self.lmp_find.grid(column=0, row=self.nlmp_array)
        self.lmp_find.insert(0, '')
        self.lmp_array_find.append(self.lmp_find)
        
        self.lmp_replace = tk.Entry(self.lmp_array_frame, width=int(1.25*self.maxwidth), font=self.font_settings)
        self.lmp_replace.grid(column=1, row=self.nlmp_array, columnspan=2)
        self.lmp_replace.insert(0, '')
        self.lmp_array_replace.append(self.lmp_replace)
        
        # adjust packing of other things in inputs frame
        self.add_button_lmp.grid(column=0, row=self.nlmp_array+1, sticky='news', columnspan=1)
        self.remove_button_lmp.grid(column=1, row=self.nlmp_array+1, sticky='news', columnspan=1)
        self.clear_button_lmp.grid(column=2, row=self.nlmp_array+1, sticky='news', columnspan=1)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.lmp_array_frame.winfo_children():
            widget.grid_configure(padx=int(xpadding/4), pady=int(ypadding/3))
        return
    
    # Function to remove last anaylsis
    def remove_last_lmp(self):
        used = []
        for index, (i, j)  in enumerate(zip(self.lmp_array_find, self.lmp_array_replace)):
            find = i.get(); replace = j.get()
            if find != '' or replace != '': used.append(index)
        if used:
            last_add = max(used)
            self.lmp_array_find[last_add].delete(0, tk.END)
            self.lmp_array_replace[last_add].delete(0, tk.END)
        else: print('No find or replace left to remove')
        return
    
    # Function to clear all files
    def clear_all_lmp(self):
        for n, i in enumerate(self.lmp_array_find):
            try: self.lmp_array_find[n].delete(0, tk.END)
            except: pass
        for n, i in enumerate(self.lmp_array_replace):
            try: self.lmp_array_replace[n].delete(0, tk.END)
            except: pass  

        self.destroy_lmp_rows(keep_first_row=True)
        return
    
    # Method to destroy lmp rows
    def destroy_lmp_rows(self, keep_first_row=False):
        """
        Destroys the dynamic peak row widgets (Entries/Comboboxes) and clears tracking lists.
        If keep_first_row=True, keeps row 1 (index 0) widgets and destroys the rest.
        """
        if keep_first_row: start = 1
        else: start = 0
    
        # destroy widgets from the end (safer)
        for idx in range(len(self.lmp_array_find)-1, start-1, -1):
            for w in (self.lmp_array_find[idx], self.lmp_array_replace[idx]):
                try: w.destroy()
                except: pass
    
        # shrink lists and update n-to-match
        self.lmp_array_find    = self.lmp_array_find[:start]
        self.lmp_array_replace = self.lmp_array_replace[:start]
        self.nlmp_array        = len(self.lmp_array_find)
        return
    
    # Function to get run mode
    def get_run_mode(self):
        self.mode = self.settings['mode']
        return
        
    # Closing command    
    def closing(self):
        print('Terminating xrd_analysis GUI'); self.root.destroy();
        return
    
    # Function to toggle if replacing log file during loading mode
    def print_selection(self):
        if (self.replace.get()) == 1:
            self.log.out(f'Will replace lmp_script when loading mode (var = {self.replace.get()})')
        else:
            self.log.out(f'Will NOT replace lmp_script when loading mode (var = {self.replace.get()})')


    # Function to get filepath for modefile
    def modefile_path(self):
        ftypes = (('Python files (.py)', '*.py'), ('all files', '*.*'))
        startpath = io_functions.path_to_string(os.path.join(self.pwd, self.modespath))
        path = filedialog.askopenfilename(initialdir=startpath, title='Open mode file?', filetypes=ftypes)
        if path:
            self.modespath = os.path.dirname(os.path.abspath(path))
            path = io_functions.path_to_string(os.path.relpath(path))
            self.modefile.delete(0, tk.END); self.modefile.insert(0, path);
            
            module = main.import_file(path)
            mode = module.mode
            self.load_mode(mode)
        else: self.log.GUI_error('ERROR could not load mode')
        
    # Function to get filepath for lmp_script
    def lmp_script_path(self):
        ftypes = (('all files', '*.*'), ('LAMMPS files (.script, .lmp)', '*.script *.lmp'))
        path = filedialog.askopenfilename(title='Open lmp_script?', filetypes=ftypes)
        if path:
            path = os.path.relpath(path)
            path = io_functions.path_to_string(path)
            self.lmp_script.delete(0, tk.END); self.lmp_script.insert(0, path)
            
            # Try getting mode from comments in lmp_script
            try:
                tmp_mode = mode_from_script.parse(path, log=None)
                if tmp_mode:
                    mode = tmp_mode
                    self.load_mode(mode)
                    lines = mode_from_script.mode2lines(mode)
                    prefix = [f'Loaded mode from: {path}', '\n\n']
                    out = prefix + lines
                    self.popup(out, title='Mode from Script', width=150)
                else:
                    print('Tried loading mode, but no mode in comments or there was a python syntax error')
            except: 
                stack_trace_string = traceback.format_exc()
                self.popup([stack_trace_string], title='Mode from Script', width=150)
                print(stack_trace_string)
        return
    
    # Function to get filepath for batch_script
    def batch_script_path(self):
        ftypes = (('Batch files (.sh)', '*.sh'), ('all files', '*.*'))
        path = filedialog.askopenfilename(title='Open batch_script?', filetypes=ftypes)
        if path:
            path = os.path.relpath(path)
            path = io_functions.path_to_string(path)
            self.batch_script.delete(0, tk.END); self.batch_script.insert(0, path)
        return
    
    # Function to get directory
    def directory_path(self):
        path = filedialog.askdirectory()
        if path:
            path = os.path.relpath(path)
            path = io_functions.path_to_string(path)
            self.parent_directory.delete(0, tk.END); self.parent_directory.insert(0, path);
        return

    # Function to save current GUI settings as a mode
    def save_mode(self):
        modes_dir = self.modespath
        if os.path.isdir(modes_dir):
            path = filedialog.asksaveasfilename(defaultextension=".py",
                                                filetypes=[("Pyton files", "*.py")],
                                                initialdir=modes_dir)
        else:
            path = filedialog.asksaveasfilename(defaultextension=".py",
                                                filetypes=[("Pyton files", "*.py")])
            
        if path:
            # Write new file and change back to pwd
            mode = self.construct_mode_dict()
            with open(path, 'w') as f:
                f.write('# -*- coding: utf-8 -*-\n')
                f.write('"""\n')
                f.write('Created by array_lmp_script.py to store a mode dictionary\n')
                f.write('called "mode", which will then allow all settings to be\n')
                f.write('loaded by clicking on this file within array_lmp_script.py.\n\n')
                #f.write('{}'.format(about.get('1.0', tk.END)))
                f.write('"""\n')
                
                f.write('# LAMMPS script find and replace. Rules:\n')
                f.write('#    find = string (recommend using $<string>)\n')
                f.write('#    replace = iterable that eval() understands\n')
                lmp_array = mode['lmp_array']
                if len(lmp_array) > 0:
                    f.write('lmp_array = {')
                    for n, find in enumerate(lmp_array):
                        string = "'{}' : '{}'".format(str(find), str(lmp_array[find]))
                        if n+1 < len(lmp_array): string += ','
                        else: string += '}'
                        if n == 0: f.write("{}\n".format(string))
                        else: f.write("{:>13}{}\n".format('', string))
                else: f.write('lmp_array = {}\n')
                f.write('\n\n')
                
                f.write('# Batch script find and replace. Rules:\n')
                f.write('#    find = string (recommend using $<string>)\n')
                f.write('#    replace = float, int, or string that eval() understands\n')
                batch_array = mode['batch_array']
                if len(batch_array) > 0:
                    f.write('batch_array = {')
                    for n, find in enumerate(batch_array):
                        string = "'{}' : '{}'".format(str(find), str(batch_array[find]))
                        if n+1 < len(batch_array): string += ','
                        else: string += '}'
                        if n == 0: f.write("{}\n".format(string))
                        else: f.write("{:>15}{}\n".format('', string))
                else: f.write('batch_array = {}\n')
                f.write('\n\n')
                
                
                f.write('# List of files to copy into the build directories (find\n')
                f.write('# and replace operations work on these file strings).\n')
                copy_files = mode['copy_files']
                if len(copy_files) > 0:
                    f.write('copy_files = [')
                    for n, file in enumerate(copy_files):
                        string = "'{}'".format(str(file))
                        if n+1 < len(copy_files): string += ','
                        else: string += ']'
                        if n == 0: f.write("{}\n".format(string))
                        else: f.write("{:>14}{}\n".format('', string))
                else: f.write('copy_files = []\n')
                f.write('\n\n')
                
                f.write('# loadable mode\n')
                f.write("mode = {")
                f.write("{}{}: '{}',\n".format('', "'lmp_script'", io_functions.path_to_string(mode['lmp_script']) ))
                f.write("{:>8}{}: '{}',\n".format('', "'batch_script'", io_functions.path_to_string(mode['batch_script']) ))
                f.write("{:>8}{}: '{}',\n".format('', "'auto_script'", mode['auto_script'] ))
                f.write("{:>8}{}: '{}',\n".format('', "'build_script'", mode['build_script'] ))
                f.write("{:>8}{}: '{}',\n".format('', "'build_directory'", mode['build_directory'] ))
                f.write("{:>8}{}: '{}',\n".format('', "'lmp_script_find'", mode['lmp_script_find'] ))
                f.write("{:>8}{}: '{}',\n".format('', "'build_script_find'", mode['build_script_find'] ))
                f.write("{:>8}{}: '{}',\n".format('', "'auto_submit'", mode['auto_submit'] ))
                f.write("{:>8}{}: {},\n".format('', "'lmp_array'", 'lmp_array' ))
                f.write("{:>8}{}: {},\n".format('', "'batch_array'", 'batch_array' ))
                f.write("{:>8}{}: {}{}\n".format('', "'copy_files'", 'copy_files', '}' ))

        self.log.out('Saving current GUI settings as a mode')
        return
    
    # Function to load mode
    def load_mode(self, mode):       
        # Start updating settings
        if self.load_replace_logfile.get() == 'True':
            self.lmp_script.delete(0, tk.END)
            self.lmp_script.insert(0, mode['lmp_script'])
            
        self.batch_script.delete(0, tk.END)
        self.batch_script.insert(0, mode['batch_script'])
            
        self.auto_script.delete(0, tk.END)
        self.auto_script.insert(0, mode['auto_script'])
        
        self.build_script.delete(0, tk.END)
        self.build_script.insert(0, mode['build_script'])
        
        self.build_directory.delete(0, tk.END)
        self.build_directory.insert(0, mode['build_directory'])
        
        self.lmp_script_find.delete(0, tk.END)
        self.lmp_script_find.insert(0, mode['lmp_script_find'])
        
        self.build_script_find.delete(0, tk.END)
        self.build_script_find.insert(0, mode['build_script_find'])
        
        self.auto_submit.delete(0, tk.END)
        self.auto_submit.insert(0, mode['auto_submit'])
        
        # Start updating lmp_array
        self.clear_all_lmp()
        lmp_array = mode['lmp_array']
        nloaded = len(self.lmp_array_find)
        for n, find in enumerate(lmp_array):
            replace = str(lmp_array[find])
            
            # add more boxes as needed
            if n >= nloaded: self.add_overloaded_analysis_lmp()
                
            self.lmp_array_find[n].delete(0, tk.END)
            self.lmp_array_find[n].insert(0, find)
            
            self.lmp_array_replace[n].delete(0, tk.END)
            self.lmp_array_replace[n].insert(0, replace)
            
        # Start updating batch_array
        self.clear_all_batch()
        batch_array = mode['batch_array']
        nloaded = len(self.batch_array_find)
        for n, find in enumerate(batch_array):
            replace = str(batch_array[find])
            
            # add more boxes as needed
            if n >= nloaded: self.add_overloaded_analysis_batch()
                
            self.batch_array_find[n].delete(0, tk.END)
            self.batch_array_find[n].insert(0, find)
            
            self.batch_array_replace[n].delete(0, tk.END)
            self.batch_array_replace[n].insert(0, replace)
            
        # Start updating batch_array
        self.clear_all_copy()
        copy_files = mode['copy_files']
        nloaded = len(self.copy_files)
        for n, file in enumerate(copy_files):
            if n >= nloaded: self.add_overloaded_analysis_copy()
            
            self.copy_files[n].delete(0, tk.END)
            self.copy_files[n].insert(0, file)
        return
    
    # Function to pop-up scrollable text
    def popup(self, out, title='Outputs', width=150):
        page = Toplevel(self.root)
        page.title(title)
        outputs = ScrolledText(page, height=30, width=width, font=('consolas', '12', 'normal'))
        outputs.pack()
        outputs.insert(tk.INSERT, '\n'.join(out))
        outputs.config(state=tk.DISABLED)
        return
    
    ########################################################################
    # Function to analyze data and plot (this is the "heart" of this code) #
    ########################################################################
    def run_main(self):
        # Generate mode dictionary from GUI options
        mode = self.construct_mode_dict()
        
        # Run main function
        main.main(mode, get_mode_from_script=False, log=self.log)
        self.popup(self.log.logged, title='Outputs')

        return