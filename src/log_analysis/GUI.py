# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.2
November 11th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.GUI_scale_settings as GUI_scale_settings
import src.log_analysis.read_log as read_log
import src.io_functions as io_functions
import src.log_analysis.main as main
from tkinter.scrolledtext import ScrolledText
from tkinter import filedialog
from tkinter import Toplevel
from tkinter import ttk
import tkinter as tk
import glob
import time
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
            mode = module.mode
        except:
            mode = {'logfile': 'UPDATE-ME',
                    'keywords': ['Step', 'Temp'],
                    'sections': 'all',
                    'xdata': '',
                    'ydata': '',
                    'xlabel': '',
                    'ylabel': '',
                    'xcompute': '',
                    'ycompute': '',
                    'analysis': [],
                    'nevery': '1',
                    'parent_directory': 'logfile'}
            self.log.GUI_error(f'ERROR loading mode file {settings["mode"]}. Internally deriving default settings. Likely launching from outside of LUNAR directory.')

        self.columns = ['Step'];
        self.mode = self.settings['mode']
        self.xdata = mode['xdata']
        self.ydata = mode['ydata']
        if self.xdata not in self.columns: self.columns.append(self.xdata)
        if self.ydata not in self.columns: self.columns.append(self.ydata)
        self.keywords = mode['keywords']
        self.sections = mode['sections']
        self.xlabel = mode['xlabel']
        self.ylabel = mode['ylabel']
        self.xcompute = mode['xcompute']
        self.ycompute = mode['ycompute']
        self.analysis = mode['analysis']
        try: self.nevery = mode['nevery']
        except: self.nevery = 1
        
        # Initialize window
        self.root = tk.Tk()
        self.root.title('LUNAR/log_analysis.py GUI v1.0')
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
        self.maxwidth = 150

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
        self.inputs_frame = tk.LabelFrame(self.frame, text='Inputs', font=self.font_settings)
        self.inputs_frame.grid(row=0, column=0, padx=self.xpadding, pady=self.ypadding)
        
        # logfile selection button
        self.logfile = tk.Entry(self.inputs_frame, width=int(1.45*self.maxwidth), font=self.font_settings)
        self.logfile.insert(0, mode['logfile'])
        self.logfile.grid(column=1, row=0, columnspan=4)
        self.logfile_button = tk.Button(self.inputs_frame, text='logfile', font=self.font_settings, command=self.logfile_path)
        self.logfile_button.grid(column=0, row=0)
                
        # logfile selection button
        self.parent_directory = tk.Entry(self.inputs_frame, width=int(1.45*self.maxwidth), font=self.font_settings)
        self.parent_directory.insert(0, mode['parent_directory'])
        self.parent_directory.grid(column=1, row=1, columnspan=4)
        self.parent_directory_button = tk.Button(self.inputs_frame, text='parent_directory', font=self.font_settings, command=self.directory_path)
        self.parent_directory_button.grid(column=0, row=1)
        
        # modes
        self.modefile = tk.Entry(self.inputs_frame, width=int(1.05*self.maxwidth), font=self.font_settings)
        self.modefile.insert(0, settings['mode'])
        self.modefile.grid(column=1, row=2)
        self.modefile_button = tk.Button(self.inputs_frame, text='mode file', font=self.font_settings, command=self.modefile_path)
        self.modefile_button.grid(column=0, row=2)        
        
        # load_replace_logfile drop down menu
        styles = [True, False]
        self.load_replace_logfile = ttk.Combobox(self.inputs_frame, values=styles, width=int(self.maxwidth/10), font=self.font_settings)
        self.load_replace_logfile.current(styles.index(settings['replace_logfile_when_loading_mode']))
        self.load_replace_logfile.grid(column=4, row=2)
        self.load_replace_logfile_label = tk.Label(self.inputs_frame, text='Replace logfile when loading mode', font=self.font_settings)
        self.load_replace_logfile_label.grid(column=3, row=2)

        
        # Add padding to all frames in self.inputs_frame
        for widget in self.inputs_frame.winfo_children():
            widget.grid_configure(padx=int(self.xpadding/2), pady=int(self.ypadding/3))
            
            
        #------------#
        # Load frame #
        #------------#
        # Initalize load data frame
        self.load_frame = tk.LabelFrame(self.frame, text='Data2load', font=self.font_settings)
        self.load_frame.grid(row=1, column=0, sticky='news', padx=self.xpadding, pady=self.ypadding)
        
        # keywords
        keywords = ','.join(self.keywords)
        self.keywords_entry = tk.Entry(self.load_frame, width=int(self.maxwidth/4.25), font=self.font_settings)
        self.keywords_entry.insert(0, keywords)
        self.keywords_entry.grid(column=0, row=1)
        self.keywords_label = tk.Label(self.load_frame, text='keywords\n(comma seperated)', font=self.font_settings)
        self.keywords_label.grid(column=0, row=0)
        
        # sections entry
        self.sections_entry = tk.Entry(self.load_frame, width=int(self.maxwidth/4), font=self.font_settings)
        self.sections_entry.insert(0, self.sections)
        self.sections_entry.grid(column=1, row=1)
        self.sections_label = tk.Label(self.load_frame, text='sections\n(all or 1 or 1,2,3 or 1-3 or 1,2,4-6 or or -1 or -2 ...)', font=self.font_settings)
        self.sections_label.grid(column=1, row=0)
        
        # xdata drop down menu
        self.xdata_dropdown = ttk.Combobox(self.load_frame, values=self.columns, width=int(self.maxwidth/6.2), font=self.font_settings,
                                              postcommand=lambda: self.xdata_dropdown.configure(values=self.columns) )
        self.xdata_dropdown.current(self.columns.index(self.xdata))
        self.xdata_dropdown.grid(column=2, row=1)
        self.xdata_dropdown_label = tk.Label(self.load_frame, text='X-data', font=self.font_settings)
        self.xdata_dropdown_label.grid(column=2, row=0)
        
        # ydata drop down menu
        self.ydata_dropdown = ttk.Combobox(self.load_frame, values=self.columns, width=int(self.maxwidth/6.2), font=self.font_settings,
                                              postcommand=lambda: self.ydata_dropdown.configure(values=self.columns) )
        self.ydata_dropdown.current(self.columns.index(self.ydata))
        self.ydata_dropdown.grid(column=3, row=1)
        self.ydata_dropdown_label = tk.Label(self.load_frame, text='Y-data', font=self.font_settings)
        self.ydata_dropdown_label.grid(column=3, row=0)
        
        # nevery
        self.nevery_entry = tk.Entry(self.load_frame, width=int(self.maxwidth/9), font=self.font_settings)
        self.nevery_entry.insert(0, self.nevery)
        self.nevery_entry.grid(column=4, row=1)
        self.nevery_label = tk.Label(self.load_frame, text='nevery\n(integer)', font=self.font_settings)
        self.nevery_label.grid(column=4, row=0)
        
        # xlabel
        self.xlabel_entry = tk.Entry(self.load_frame, width=int(self.maxwidth/4.25), font=self.font_settings)
        self.xlabel_entry.insert(0, self.xlabel)
        self.xlabel_entry.grid(column=5, row=1)
        self.xlabel_label = tk.Label(self.load_frame, text='X-label\n(string)', font=self.font_settings)
        self.xlabel_label.grid(column=5, row=0)
        
        # ylabel
        self.ylabel_entry = tk.Entry(self.load_frame, width=int(self.maxwidth/4.25), font=self.font_settings)
        self.ylabel_entry.insert(0, self.ylabel)
        self.ylabel_entry.grid(column=6, row=1)
        self.ylabel_label = tk.Label(self.load_frame, text='Y-label\n(string)', font=self.font_settings)
        self.ylabel_label.grid(column=6, row=0)
        
        # Add padding to all frames in self.load_frame
        for widget in self.load_frame.winfo_children():
            widget.grid_configure(padx=int(self.xpadding/2), pady=int(self.ypadding/3))
            
        #----------------#
        # Computes frame #
        #----------------#
        # Initalize computes frame
        self.computes_frame = tk.LabelFrame(self.frame, text='Computes', font=self.font_settings)
        self.computes_frame.grid(row=2, column=0, sticky='news', padx=self.xpadding, pady=self.ypadding)
            
        # xcompute
        self.xcompute_entry = tk.Entry(self.computes_frame, width=int(1.3*self.maxwidth), font=self.font_settings)
        self.xcompute_entry.insert(0, self.xcompute)
        self.xcompute_entry.grid(column=1, row=1)
        self.xcompute_label = tk.Label(self.computes_frame, text='X-compute', font=self.font_settings)
        self.xcompute_label.grid(column=0, row=1)
        
        # ycompute
        self.ycompute_entry = tk.Entry(self.computes_frame, width=int(1.3*self.maxwidth), font=self.font_settings)
        self.ycompute_entry.insert(0, self.ycompute)
        self.ycompute_entry.grid(column=1, row=2)
        self.ycompute_label = tk.Label(self.computes_frame, text='Y-compute', font=self.font_settings)
        self.ycompute_label.grid(column=0, row=2)
        
        # Button to load compute help
        self.compute_help = tk.Button(self.computes_frame, text='compute\nhelp', font=self.font_settings, width=int(self.maxwidth/8), command=self.compute_options)
        self.compute_help.grid(column=3, row=1, columnspan=1, rowspan=2)
        
        # Add padding to all frames in self.computes_frame
        for widget in self.computes_frame.winfo_children():
            widget.grid_configure(padx=int(self.xpadding/2), pady=int(self.ypadding/3))
            
        #----------------#
        # Analysis frame #
        #----------------#
        # Initalize analysis frame
        self.analysis_frame = tk.LabelFrame(self.frame, text='Analysis options', font=self.font_settings)
        self.analysis_frame.grid(row=3, column=0, sticky='news', padx=self.xpadding, pady=self.ypadding)
        
        # file selection labels
        self.method_label = tk.Label(self.analysis_frame, text='method', font=self.font_settings)
        self.method_label.grid(column=0, row=0)
        self.xlo_label = tk.Label(self.analysis_frame, text='X-lo (float or empty)', font=self.font_settings)
        self.xlo_label.grid(column=1, row=0)
        self.xhi_label = tk.Label(self.analysis_frame, text='X-hi (float or empty)', font=self.font_settings)
        self.xhi_label.grid(column=2, row=0)
        self.misc_label = tk.Label(self.analysis_frame, text='Misc settings (string or empty)', font=self.font_settings)
        self.misc_label.grid(column=3, row=0)
        self.name_label = tk.Label(self.analysis_frame, text='Name (string or empty)', font=self.font_settings)
        self.name_label.grid(column=4, row=0)

        
        # file selection button and qty
        self.nanalysis = len(self.analysis)
        self.supported_methods = ['average', 'linear regression', 'moving average', 'hyperbola', 'piecewise-regression', 'cursor', 'skip',
                                  'spline-integration', 'Whittaker-Eilers', 'minimum', 'maximum', 'Butterworth (low pass)', 'iFFT filter', 'LOWESS',
                                  'Regression Fringe Response Modulus', 'LAMMPS data (remove from plot)', 'LAMMPS data (apply moving average)',
                                  'LAMMPS data (apply Butterworth filter)', 'LAMMPS data (apply Whittaker-Eilers)', 'LAMMPS data (fit polynomial)',
                                  'LAMMPS data (LOWESS)', 'LAMMPS data (X-sort)', 'LAMMPS data (apply iFFT filter)', 'write plotted data to csv file',
                                  'Calculus: Differentiate Data', 'Calculus: Integrate Data', 'Regression Fringe Response Thermal']
        self.supported_methods = sorted(self.supported_methods, key=lambda x: x[0].lower()) # sort list by first letter of each method (x[0].lower())
        self.methods = []; self.xlos = []; self.xhis = []; self.miscs = []; self.names = [];
        for n in range(1, self.nanalysis+1):
            method, xlo, xhi, misc, name = self.analysis[n-1]
            if method == '': break
            if name == '': name = 'analysis-{}'.format(n)
            
            if method not in self.supported_methods:
                self.supported_methods.append(method)
            
            self.method = ttk.Combobox(self.analysis_frame, values=self.supported_methods, width=int(self.maxwidth/4), font=self.font_settings)
            self.method.current(self.supported_methods.index(method))
            self.method.grid(column=0, row=n)
            self.methods.append(self.method)
            
            self.xlo = tk.Entry(self.analysis_frame, width=int(self.maxwidth/6), font=self.font_settings)
            self.xlo.grid(column=1, row=n)
            self.xlo.insert(0, xlo)
            self.xlos.append(self.xlo)
            
            self.xhi = tk.Entry(self.analysis_frame, width=int(self.maxwidth/6), font=self.font_settings)
            self.xhi.grid(column=2, row=n)
            self.xhi.insert(0, xhi)
            self.xhis.append(self.xhi)
            
            self.misc = tk.Entry(self.analysis_frame, width=int(self.maxwidth/1.5), font=self.font_settings)
            self.misc.grid(column=3, row=n)
            self.misc.insert(0, misc)
            self.miscs.append(self.misc)
            
            self.name = tk.Entry(self.analysis_frame, width=int(self.maxwidth/4), font=self.font_settings)
            self.name.grid(column=4, row=n)
            self.name.insert(0, name)
            self.names.append(self.name)
            
        # Button to add a file
        self.add_button = tk.Button(self.analysis_frame, text='add analysis to stack', font=self.font_settings, width=int(self.maxwidth/4.5), command=self.add2stack)
        self.add_button.grid(column=0, row=self.nanalysis+1, columnspan=1)
            
        # Button to remove a file
        self.remove_button = tk.Button(self.analysis_frame, text='remove last analysis from stack', font=self.font_settings, width=int(self.maxwidth/3.5), command=self.remove_last)
        self.remove_button.grid(column=1, row=self.nanalysis+1, sticky='news', columnspan=2)
        
        # Button to clear all files
        self.clear_button = tk.Button(self.analysis_frame, text='clear stack', font=self.font_settings, width=int(self.maxwidth/2.1), command=self.clear_all)
        self.clear_button.grid(column=3, row=self.nanalysis+1, columnspan=1)
        
        # Button load compute help
        self.analysis_help = tk.Button(self.analysis_frame, text='analysis help', font=self.font_settings, width=int(self.maxwidth/5), command=self.analysis_options)
        self.analysis_help.grid(column=4, row=self.nanalysis+1, columnspan=1)
            
        # Add padding to all frames in self.analysis_frame
        for widget in self.analysis_frame.winfo_children():
            widget.grid_configure(padx=int(self.xpadding/4), pady=int(self.ypadding/3))
            
        
        #-------------------------#
        # Initalize execute frame #
        #-------------------------#
        self.execute_frame = tk.LabelFrame(self.frame, borderwidth=0, highlightthickness=0)
        self.execute_frame.grid(column=0, row=4, padx=self.xpadding, pady=self.ypadding)
                
        # data2csv
        self.data2csv_btn = tk.Button(self.execute_frame, width=int(self.maxwidth/2.25), text='write loaded data to csv', command=self.write2csv, font=self.font_settings)
        self.data2csv_btn.grid(column=0, row=0)
        
        # save_mode
        self.save_mode_btn = tk.Button(self.execute_frame, width=int(self.maxwidth/2.25), text='save settings as mode', command=self.save_mode, font=self.font_settings)
        self.save_mode_btn.grid(column=1, row=0)
        
        # update plot
        self.update_plot_btn = tk.Button(self.execute_frame, width=int(self.maxwidth/2.25), text='update plot', command=self.analyze_and_plot, font=self.font_settings)
        self.update_plot_btn.grid(column=2, row=0)
        
        # Add padding to all frames in self.analysis_frame
        for widget in self.execute_frame.winfo_children():
            widget.grid_configure(padx=int(self.xpadding/2), pady=int(self.ypadding/3))
        
        #------------------------#
        # Run mainloop and close #
        #------------------------#
        self.root.protocol('WM_DELETE_WINDOW', self.closing)
        self.root.mainloop()
    
    ################################
    # Functions to call as methods #
    ################################
    # Method to get basename of file and build directories if needed
    def get_basename_and_builder_dirs(self):
        # Get currently defined parent_directory variable
        logfile = self.logfile.get()
        parent_directory = self.parent_directory.get()
        
        # Find present working directory
        pwd = os.getcwd()
        
        # Find/create paths to store code results
        path = os.path.join(pwd, parent_directory)
        
        # Going to use io_functions get_dir_from_topofile() function which uses 'topofile'
        # string not 'logfile' string. So convert 'logfile' to 'topofile' if applicable
        if 'logfile' in parent_directory:
            parent_directory = parent_directory.replace('logfile', 'topofile')
            self.log.out('Using path from logfile to set parent_directory ...')
            path = io_functions.get_dir_from_topofile(logfile, parent_directory)
            
        # Check if path exists. IF not create.
        if not os.path.isdir(path):
            os.makedirs(path, exist_ok=True)
            
        # Set basename 
        root = os.path.basename(logfile)
        basename = os.path.join(path, root)
        return basename
    
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
        self.nanalysis += 1
        self.method = ttk.Combobox(self.analysis_frame, values=self.supported_methods, width=int(self.maxwidth/4), font=self.font_settings)
        self.method.current(self.supported_methods.index('skip'))
        self.method.grid(column=0, row=self.nanalysis)
        self.methods.append(self.method)
        
        self.xlo = tk.Entry(self.analysis_frame, width=int(self.maxwidth/6), font=self.font_settings)
        self.xlo.grid(column=1, row=self.nanalysis)
        #self.xlo.insert(0, 0)
        self.xlos.append(self.xlo)
        
        self.xhi = tk.Entry(self.analysis_frame, width=int(self.maxwidth/6), font=self.font_settings)
        self.xhi.grid(column=2, row=self.nanalysis)
        #self.xhi.insert(0, 1)
        self.xhis.append(self.xhi)
        
        self.misc = tk.Entry(self.analysis_frame, width=int(self.maxwidth/1.5), font=self.font_settings)
        self.misc.grid(column=3, row=self.nanalysis)
        self.misc.insert(0, '')
        self.miscs.append(self.misc)
        
        self.name = tk.Entry(self.analysis_frame, width=int(self.maxwidth/4), font=self.font_settings)
        self.name.grid(column=4, row=self.nanalysis)
        self.name.insert(0, '')
        self.names.append(self.name)
        
        # adjust packing of other things in inputs frame
        self.add_button.grid(column=0, row=self.nanalysis+1, columnspan=1)
        self.remove_button.grid(column=1, row=self.nanalysis+1, sticky='news', columnspan=2)
        self.clear_button.grid(column=3, row=self.nanalysis+1, columnspan=1)
        self.analysis_help.grid(column=4, row=self.nanalysis+1, columnspan=1)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.analysis_frame.winfo_children():
            widget.grid_configure(padx=int(xpadding/4), pady=int(ypadding/3))
        return
    
    # Function to remove last anaylsis
    def remove_last(self):
        used = []
        for n, (i, j, k, l, m)  in enumerate(zip(self.methods, self.xlos, self.xhis, self.miscs, self.names)):
            method = i.get(); xlo = j.get(); xhi = k.get(); misc = l.get(); name = m.get();
            if method != '' or xlo != '' or xhi != '' or name != '': used.append(n)
        if used:
            last_add = max(used)
            self.methods[last_add].current(self.supported_methods.index('skip'))
            self.xlos[last_add].delete(0, tk.END)
            self.xhis[last_add].delete(0, tk.END)
            self.miscs[last_add].delete(0, tk.END)
            self.names[last_add].delete(0, tk.END)
        else: print('No files or tags left to remove')
        return
    
    # Function to clear all files
    def clear_all(self):
        for n, i in enumerate(self.methods):
            try: self.methods[n].current(self.supported_methods.index('skip'))
            except: pass
        for n, i in enumerate(self.xlos):
            try: self.xlos[n].delete(0, tk.END)
            except: pass
        for n, i in enumerate(self.xhis):
            try: self.xhis[n].delete(0, tk.END)
            except: pass
        for n, i in enumerate(self.miscs):
            try: self.miscs[n].delete(0, tk.END)
            except: pass
        for n, i in enumerate(self.names):
            try: self.names[n].delete(0, tk.END)
            except: pass
        return
    
    # Function to get filepath for logfile
    def logfile_path(self):
        ftypes = (('all files', '*.*'), ('LAMMPS log files (.lammps, .log, .txt)', '*.lammps *.log *.txt'))
        path = filedialog.askopenfilename(title='Open logfile?', filetypes=ftypes)
        if path:
            self.filepath = os.path.dirname(os.path.abspath(path))
            path = os.path.relpath(path)
            self.logfile.delete(0, tk.END); self.logfile.insert(0, path)
            
            # Update self.columns
            if os.path.isfile(path):
                try: keywords = self.keywords.get().split(',')
                except: keywords = ['Step', 'Temp', 'Density', 'Press', 'TotEng', 'KinEng', 'PotEng', 'Lx', 'Ly', 'Lz']
                try: 
                    sections = self.sections_entry.get()
                    if sections == '': sections = 'all'
                except: sections = 'all'
                log = read_log.file(path, keywords=keywords)
                data = log.get_data(sections, remove_duplicates=True, pflag=True) # {column-name:[lst of data]}
                self.columns = list(data.keys())
        return
    
    # Function to get directory
    def directory_path(self):
        path = filedialog.askdirectory()
        if path:
            path = os.path.relpath(path)
            self.parent_directory.delete(0, tk.END); self.parent_directory.insert(0, path);
        return
    
    # Function to write data to csv
    def write2csv(self):
        # Read logfile and get data
        try:
            logfile = self.logfile.get()
            if os.path.isfile(logfile):
                keywords = self.keywords_entry.get().split(',')
                sections = self.sections_entry.get()
                log = read_log.file(logfile, keywords=keywords)
                data2write = log.get_data(sections, remove_duplicates=True, pflag=True) # {column-name:[lst of data]}
                self.columns = list(data2write.keys()) # update columns to push to xdata, ydata drop downs
            else: data2write = {}; print(f'lammps logfile {logfile} does not exist');
        except: data2write = {}; print('ERROR failed to read logfile and load data')
        
        # writedata to csv
        if data2write:
            try: nevery = int(self.nevery_entry.get())
            except: nevery = 1
            if nevery > 1:
                for column in data2write:
                    tmp = data2write[column][::nevery]
                    data2write[column] = tmp
            
            csvname = '{}.csv'.format(self.get_basename_and_builder_dirs())
            print(f'Writing {csvname}')
            with open(csvname, 'w') as f:
                # invert data
                ncolumns = len(data2write); nrows = max(map(len, list(data2write.values())))
                matrix = [[0]*ncolumns for n in range(nrows)]
                titles = sorted(data2write.keys())
                for i in range(nrows):
                    for j, name in enumerate(titles):
                        #print(i, j, name, data2write[name][i])
                        matrix[i][j] = str(data2write[name][i])
                
                # Join with comma's and write titles
                titles = ', '.join(titles);
                f.write('{}\n'.format(titles));
                
                # write rows
                for row in matrix:
                    row = ', '.join(row);
                    f.write('{}\n'.format(row));
        else: print('ERROR could not write data to csv file')
        return
    
    # Function to get run mode
    def get_run_mode(self):
        self.mode = self.settings['mode']
        return
        
    # Closing command    
    def closing(self):
        print('Terminating log_analysis GUI'); self.root.destroy();
        return
    
    # Function to toggle if replacing log file during loading mode
    def print_selection(self):
        if (self.replace.get()) == 1:
            self.log.out(f'Will replace logfile when loading mode (var = {self.replace.get()})')
        else:
            self.log.out(f'Will NOT replace logfile when loading mode (var = {self.replace.get()})')
    
    # Analysis quick help page
    def analysis_options(self):
        try: # Try to get text from GUI_help_page.txt file
            txt = os.path.join(self.pwd, 'src/GUI_quick_help_pages/log_analysis_methods.txt')
            logged = []
            with open(txt, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    # Strip comment's and split by whitespace
                    line = line.split('#')[0]
                    line = line.rstrip()
                    logged.append(line)
        except: # except something failed
            logged.append('FAILED to read LUNAR/src/GUI_quick_help_pages/log_analysis_methods..txt document.')
            logged.append('Most likely cause is the log_analysis_methods.txt file was renamed in the')
            logged.append('LUNAR/src/GUI_quick_help_pages/ directory or directory names were changed.')
        self.popup(logged, title='Analysis help')
        return

    # Compute quick help page
    def compute_options(self):
        try: # Try to get text from GUI_help_page.txt file
            txt = os.path.join(self.pwd, 'src/GUI_quick_help_pages/log_analysis_computes.txt')
            logged = []
            with open(txt, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    # Strip comment's and split by whitespace
                    line = line.split('#')[0]
                    line = line.rstrip()
                    logged.append(line)
        except: # except something failed
            logged.append('FAILED to read LUNAR/src/GUI_quick_help_pages/log_analysis_computes.txt document.')
            logged.append('Most likely cause is the log_analysis_computes.txt file was renamed in the')
            logged.append('LUNAR/src/GUI_quick_help_pages/ directory or directory names were changed.')
        self.popup(logged, title='Compute help')
        return
            
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
        
    # Function to save current GUI settings as a mode
    def save_mode(self):
        # Function to get directory
        def directory_path():
            path = filedialog.askdirectory(initialdir=self.pwd)
            if path:
                path = os.path.relpath(path)
                parent_directory.delete(0, tk.END); parent_directory.insert(0, path);
        
        # save_mode print(__doc__)
        def save():
            if filename.get() == '':
                self.log.GUI_error('ERROR attempting to save GUI settings as mode, but filename is blank.')
            elif filename.get().count(' ') > 0:
                self.log.GUI_error(f'ERROR attempting to save GUI settings as mode, but filename "{filename.get()}" has whitespaces. Use "_" to seperate words.')
            else:
                # setup path and where to write
                path = os.path.join(self.pwd, parent_directory.get())
                name = '{}.py'.format(filename.get())
                self.log.out(f'"{name}" will be saved at {path}')
                if not os.path.isdir(path):
                    os.makedirs(path, exist_ok=True)
                os.chdir(path)
                
                # Write new file and change back to pwd
                with open(name, 'w') as f:
                    f.write('# -*- coding: utf-8 -*-\n')
                    f.write('"""\n')
                    f.write('Created by log_analysis.py to store a mode dictionary\n')
                    f.write('called "mode", which will then allow all settings to be\n')
                    f.write('loaded by clicking on this file within log_anaylsis.py.\n\n')
                    f.write('{}'.format(about.get('1.0', tk.END)))
                    f.write('"""\n')
                    
                    f.write('# analysis list\n')
                    if len(self.methods) > 0:
                        analysis = [] # [method, xlo, xhi, miscs, names]
                        for i, j, k, l, m  in zip(self.methods, self.xlos, self.xhis, self.miscs, self.names):
                            try: 
                                method = i.get(); xlo = j.get(); xhi = k.get(); misc = l.get(); name = m.get();
                                if method != 'skip':
                                    analysis.append([method, xlo, xhi, misc, name])
                            except: pass
                        if analysis:
                            f.write('analysis = [')
                            for n, (method, xlo, xhi, misc, name) in enumerate(analysis, 1):
                                if method == 'skip': continue
                                try:
                                    if xlo == '': xlo = "'{}'".format(xlo)
                                    if xhi == '': xhi = "'{}'".format(xhi)
                                    string = "['{}', {}, {}, '{}', '{}']".format(method, xlo, xhi, misc, name)
                                    
                                    if n < len(analysis): string += ','
                                    else: string += ']'
                                    
                                    if n == 1: f.write("{}\n".format(string))
                                    else: f.write("{:>12}{}\n".format('', string))
                                except: pass
                        else: f.write('analysis = []\n')
                    else: f.write('analysis = []\n')
                    f.write('\n')
                    
                    f.write('# loadable mode\n')
                    f.write("mode = {")
                    f.write("{}{}: '{}',\n".format('', "'logfile'", io_functions.path_to_string(self.logfile.get())))
                    keywords = ["'{}'".format(str(i)) for i in self.keywords_entry.get().split(',')]
                    f.write("{:>8}{}: [{}],\n".format('',"'keywords'", ', '.join(keywords)))
                    f.write("{:>8}{}: '{}',\n".format('',"'sections'", self.sections_entry.get()))
                    
                    if "'" in self.xdata_dropdown.get():
                        f.write('{:>8}{}: "{}",\n'.format('', '"xdata"', self.xdata_dropdown.get()))
                    else:
                        f.write("{:>8}{}: '{}',\n".format('',"'xdata'", self.xdata_dropdown.get()))
                    
                    if "'" in self.ydata_dropdown.get():
                        f.write('{:>8}{}: "{}",\n'.format('', '"ydata"', self.ydata_dropdown.get()))
                    else:    
                        f.write("{:>8}{}: '{}',\n".format('',"'ydata'", self.ydata_dropdown.get()))
                    
                    f.write("{:>8}{}: '{}',\n".format('',"'xlabel'", self.xlabel_entry.get()))
                    f.write("{:>8}{}: '{}',\n".format('',"'ylabel'", self.ylabel_entry.get()))
                    f.write("{:>8}{}: '{}',\n".format('',"'xcompute'", self.xcompute_entry.get()))
                    f.write("{:>8}{}: '{}',\n".format('',"'ycompute'", self.ycompute_entry.get()))
                    f.write("{:>8}{}: {},\n".format('',"'analysis'", 'analysis'))
                    f.write("{:>8}{}: '{}',\n".format('',"'nevery'", self.nevery_entry.get()))
                    f.write("{:>8}{}: '{}'\n".format('',"'parent_directory'", io_functions.path_to_string(self.parent_directory.get())))
                    f.write("{:>8}{}\n\n".format('', '}'))
                os.chdir(self.pwd)
            return
        
        self.log.out('Saving current GUI settings as a mode')
        save_frame = Toplevel(self.root)
        save_frame.title('Save GUI settings as mode')
        save_frame.resizable(width=False, height=False)
        
        filename = tk.Entry(save_frame, width=67, font=self.font_settings)
        filename.insert(0, 'FML_analysis_description')
        filename.bind("<Button-1>", lambda a: filename.delete(0, tk.END))
        filename.grid(column=1, row=0)
        filename_label = tk.Label(save_frame, text='filename (no .py extension)', font=self.font_settings)
        filename_label.grid(column=0, row=0)
        
        about = tk.Text(save_frame, width=50, height=4)
        about.insert(tk.END, 'Date: XX/YY/ZZZZ\n')
        about.insert(tk.END, 'Author: John Doe\n')
        about.insert(tk.END, 'Purpose: Mode to find xxx\n')
        about.grid(column=1, row=1)
        about_label = tk.Label(save_frame, text='about mode', font=self.font_settings)
        about_label.grid(column=0, row=1)
        
        parent_directory = tk.Entry(save_frame, width=67, font=self.font_settings)
        parent_directory.insert(0, self.modespath)
        parent_directory.grid(column=1, row=2)
        dir_button = tk.Button(save_frame, text='parent_directory', font=self.font_settings, command=directory_path)
        dir_button.grid(column=0, row=2)
        
        save_btn = tk.Button(save_frame, width=85, text='save as mode', command=save, font=self.font_settings)
        save_btn.grid(column=0, row=3, columnspan=3)        

        # Add padding to all frames in self.analysis_frame
        for widget in save_frame.winfo_children():
            widget.grid_configure(padx=self.xpadding, pady=int(self.ypadding/3))
        save_frame.mainloop()
    
    # Function to load mode
    def load_mode(self, mode):       
        # Start updating settings
        if self.load_replace_logfile.get() == 'True':
            self.logfile.delete(0, tk.END)
            self.logfile.insert(0, mode['logfile'])
            
        self.parent_directory.delete(0, tk.END)
        self.parent_directory.insert(0, mode['parent_directory'])
        
        self.keywords = ','.join(mode['keywords'])
        self.keywords_entry.delete(0, tk.END)
        self.keywords_entry.insert(0, self.keywords)
    
        self.sections = mode['sections']        
        self.sections_entry.delete(0, tk.END)
        self.sections_entry.insert(0, self.sections)

        self.xdata = mode['xdata']
        self.ydata = mode['ydata']
        if self.xdata not in self.columns: self.columns.append(self.xdata)
        if self.ydata not in self.columns: self.columns.append(self.ydata)
        self.xdata_dropdown.config(values=self.columns)
        self.ydata_dropdown.config(values=self.columns)
        self.xdata_dropdown.current(self.columns.index(self.xdata))
        self.ydata_dropdown.current(self.columns.index(self.ydata))
        
        self.xcompute = mode['xcompute']
        self.ycompute = mode['ycompute']
        self.xcompute_entry.delete(0, tk.END)
        self.ycompute_entry.delete(0, tk.END)
        self.xcompute_entry.insert(0, self.xcompute)
        self.ycompute_entry.insert(0, self.ycompute)
        
        self.xlabel = mode['xlabel']
        self.ylabel = mode['ylabel']
        self.xlabel_entry.delete(0, tk.END)
        self.ylabel_entry.delete(0, tk.END)
        self.xlabel_entry.insert(0, self.xlabel)
        self.ylabel_entry.insert(0, self.ylabel)
        
        try: self.nevery = mode['nevery']
        except: self.nevery = '1'
        self.nevery_entry.delete(0, tk.END)
        self.nevery_entry.insert(0, self.nevery)
        
        # Start updating analysis
        self.clear_all()
        nloaded = len(self.methods)
        analysis = mode['analysis']
        for n, i in enumerate(analysis):
            method, xlo, xhi, misc, name = i
            if n >= nloaded: # add more boxes as needed
                self.add_overloaded_analysis()
            
            if method not in self.supported_methods:
                self.supported_methods.append(method)
            try: 
                index = self.supported_methods.index(method)
                self.methods[n].current(index)
            except: 
                self.log.GUI_error(f'ERROR could not find proper index of {method}. Did not update method')
                
            self.xlos[n].delete(0, tk.END)
            self.xlos[n].insert(0, xlo)
            self.xhis[n].delete(0, tk.END)
            self.xhis[n].insert(0, xhi)
            self.miscs[n].delete(0, tk.END)
            self.miscs[n].insert(0, misc)
            self.names[n].delete(0, tk.END)
            self.names[n].insert(0, name)
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
    def analyze_and_plot(self):
        # Generate mode dictionary from GUI options
        mode = {}
        mode['logfile'] = self.logfile.get()
        mode['keywords'] = self.keywords_entry.get().split(',')
        mode['sections'] = self.sections_entry.get()
        mode['xdata'] = self.xdata_dropdown.get()
        mode['ydata'] = self.ydata_dropdown.get()
        mode['xlabel'] = self.xlabel_entry.get()
        mode['ylabel'] = self.ylabel_entry.get()
        mode['xcompute'] = self.xcompute_entry.get()
        mode['ycompute'] = self.ycompute_entry.get()
        mode['nevery'] = self.nevery_entry.get()
        mode['parent_directory'] = self.parent_directory.get()
        mode['analysis'] = []
        for n, (i, j, k, l, m)  in enumerate(zip(self.methods, self.xlos, self.xhis, self.miscs, self.names)):
            analysis = [i.get(), j.get(), k.get(), l.get(), m.get()]
            mode['analysis'].append(analysis)
        
        # Analyze results with main.analysis class
        dpi = self.settings['image-dpi']
        savefig = self.settings['save-fig']
        
        # Set up Tristan's "array" analysis using recursion
        if not os.path.isfile(str(mode['logfile'])):
            files = glob.glob(mode['logfile']); array_time = time.time(); analyzed = None; 
            outputs = {} # {filename : dict-of-outputs}
            if files:
                for n, file in enumerate(files, 1):
                    self.log.clear_all()
                    self.log.out('\n\nUsing array input option:')
                    self.log.out(' - logfile      : {}'.format(mode['logfile']))
                    self.log.out(' - matched file : {}'.format(file))
                    self.log.out(' - elapsed time : {:.2f} (seconds)'.format(time.time() - array_time))
                    self.log.out(' - progress     : {} of {} ({:.2f}%)'.format(n, len(files), 100*(n/len(files))))
                    if file.endswith('.jpeg') or file.endswith('.log.lunar') or file.endswith('.eps'):
                        self.log.warn(f' - WARNING matched file {file} has an extension that does not make sense to process. Skipping file')
                        continue
                    mode['logfile'] = file
                    try: # we dont want crashes to exit this loop
                        analyzed = main.analysis(mode, plot=True, savefig=savefig, dpi=dpi, log=self.log, log_clear=False)
                        outputs[file] = analyzed.outputs
                    except: pass
            print('\a') # Alert
            # if outputs:
            #     self.array_csv(outputs)
        else:
            analyzed = main.analysis(mode, plot=True, savefig=savefig, dpi=dpi, log=self.log)
            self.columns = analyzed.columns # Update columns
            #self.popup(self.log.logged, title='Outputs', width=150)
        return
    
    def array_csv(self, outputs):
        data2skip = set(['Raw-data'])
        for file in outputs:
            print()
            print(file)
            for data in outputs[file]:
                if data in data2skip: continue
                if len(outputs[file][data]) > 10: continue
                print(data, outputs[file][data])
        return