# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
August 13th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.log_analysis.misc_funcs as misc_funcs
import src.log_analysis.read_log as read_log
import src.io_functions as io_functions
from tkinter.scrolledtext import ScrolledText
import matplotlib.pyplot as plt
from tkinter import filedialog
from tkinter import Toplevel
from tkinter import ttk
import tkinter as tk
import numpy as np
import math
import time
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
        self.filepath = self.pwd
        self.modespath = settings['modes-dir']
        
        # Set defaults
        self.settings = settings
        module = misc_funcs.import_file(settings['mode'])
        mode = module.mode
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
        self.xpadding = 20
        self.ypadding = 10
        self.maxwidth = 140
        
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
        
        # modes
        self.modefile = tk.Entry(self.inputs_frame, width=self.maxwidth, font=self.font_settings)
        self.modefile.insert(0, settings['mode'])
        self.modefile.grid(column=1, row=1)
        self.modefile_button = tk.Button(self.inputs_frame, text='mode file', font=self.font_settings, command=self.modefile_path)
        self.modefile_button.grid(column=0, row=1)        
        
        # load_replace_logfile drop down menu
        styles = [True, False]
        self.load_replace_logfile = ttk.Combobox(self.inputs_frame, values=styles, width=int(self.maxwidth/10), font=self.font_settings)
        self.load_replace_logfile.current(styles.index(settings['replace_logfile_when_loading_mode']))
        self.load_replace_logfile.grid(column=4, row=1)
        self.load_replace_logfile_label = tk.Label(self.inputs_frame, text='Replace logfile when loading mode', font=self.font_settings)
        self.load_replace_logfile_label.grid(column=3, row=1)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.inputs_frame.winfo_children():
            widget.grid_configure(padx=self.xpadding, pady=int(self.ypadding/3))
            
            
        #------------#
        # Load frame #
        #------------#
        # Initalize load data frame
        self.load_frame = tk.LabelFrame(self.frame, text='Data2load', font=self.font_settings)
        self.load_frame.grid(row=1, column=0, sticky='news', padx=self.xpadding, pady=self.ypadding)
        
        # keywords
        keywords = ','.join(self.keywords)
        self.keywords_entry = tk.Entry(self.load_frame, width=int(self.maxwidth/5), font=self.font_settings)
        self.keywords_entry.insert(0, keywords)
        self.keywords_entry.grid(column=0, row=1)
        self.keywords_label = tk.Label(self.load_frame, text='keywords\n(comma seperated)', font=self.font_settings)
        self.keywords_label.grid(column=0, row=0)
        
        # sections entry
        self.sections_entry = tk.Entry(self.load_frame, width=int(self.maxwidth/3), font=self.font_settings)
        self.sections_entry.insert(0, self.sections)
        self.sections_entry.grid(column=1, row=1)
        self.sections_label = tk.Label(self.load_frame, text='sections\n(all or 1 or 1,2,3 or 1-3 or 1,2,4-6 or ...)', font=self.font_settings)
        self.sections_label.grid(column=1, row=0)
        
        # xdata drop down menu
        self.xdata_dropdown = ttk.Combobox(self.load_frame, values=self.columns, width=int(self.maxwidth/6), font=self.font_settings,
                                              postcommand=lambda: self.xdata_dropdown.configure(values=self.columns) )
        self.xdata_dropdown.current(self.columns.index(self.xdata))
        self.xdata_dropdown.grid(column=2, row=1)
        self.xdata_dropdown_label = tk.Label(self.load_frame, text='X-data', font=self.font_settings)
        self.xdata_dropdown_label.grid(column=2, row=0)
        
        # ydata drop down menu
        self.ydata_dropdown = ttk.Combobox(self.load_frame, values=self.columns, width=int(self.maxwidth/6), font=self.font_settings,
                                              postcommand=lambda: self.ydata_dropdown.configure(values=self.columns) )
        self.ydata_dropdown.current(self.columns.index(self.ydata))
        self.ydata_dropdown.grid(column=3, row=1)
        self.ydata_dropdown_label = tk.Label(self.load_frame, text='Y-data', font=self.font_settings)
        self.ydata_dropdown_label.grid(column=3, row=0)
        
        # xlabel
        self.xlabel_entry = tk.Entry(self.load_frame, width=int(self.maxwidth/4), font=self.font_settings)
        self.xlabel_entry.insert(0, self.xlabel)
        self.xlabel_entry.grid(column=4, row=1)
        self.xlabel_label = tk.Label(self.load_frame, text='X-label\n(string)', font=self.font_settings)
        self.xlabel_label.grid(column=4, row=0)
        
        # ylabel
        self.ylabel_entry = tk.Entry(self.load_frame, width=int(self.maxwidth/4), font=self.font_settings)
        self.ylabel_entry.insert(0, self.ylabel)
        self.ylabel_entry.grid(column=5, row=1)
        self.ylabel_label = tk.Label(self.load_frame, text='Y-label\n(string)', font=self.font_settings)
        self.ylabel_label.grid(column=5, row=0)
        
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
        self.compute_help = tk.Button(self.computes_frame, text='compute\nhelp', font=self.font_settings, width=int(self.maxwidth/9), command=self.compute_options)
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
        self.supported_methods = ['average', 'linear regression', 'moving average', 'hyperbola', 'piecewise-regression', 
                                  'spline-integration', 'cursor', 'minimum', 'maximum', 'butterworth (low pass)',
                                  'LAMMPS data (remove from plot)',
                                  'LAMMPS data (apply moving average)', 'LAMMPS data (apply butterworth filter)', 'skip']
        self.supported_methods = sorted(self.supported_methods, key=lambda x: x[0].lower()) # sort list by first letter of each method (x[0].lower())
        self.methods = []; self.xlos = []; self.xhis = []; self.miscs = []; self.names = [];
        for n in range(1, self.nanalysis+1):
            method, xlo, xhi, misc, name = self.analysis[n-1]
            if method == '': break
            if name == '': name = 'analysis-{}'.format(n)
            
            self.method = ttk.Combobox(self.analysis_frame, values=self.supported_methods, width=int(self.maxwidth/4), font=self.font_settings)
            self.method.current(self.supported_methods.index(method))
            self.method.grid(column=0, row=n)
            self.methods.append(self.method)
            
            self.xlo = tk.Entry(self.analysis_frame, width=int(self.maxwidth/5), font=self.font_settings)
            self.xlo.grid(column=1, row=n)
            self.xlo.insert(0, xlo)
            self.xlos.append(self.xlo)
            
            self.xhi = tk.Entry(self.analysis_frame, width=int(self.maxwidth/5), font=self.font_settings)
            self.xhi.grid(column=2, row=n)
            self.xhi.insert(0, xhi)
            self.xhis.append(self.xhi)
            
            self.misc = tk.Entry(self.analysis_frame, width=int(self.maxwidth/2), font=self.font_settings)
            self.misc.grid(column=3, row=n)
            self.misc.insert(0, misc)
            self.miscs.append(self.misc)
            
            self.name = tk.Entry(self.analysis_frame, width=int(self.maxwidth/3), font=self.font_settings)
            self.name.grid(column=4, row=n)
            self.name.insert(0, name)
            self.names.append(self.name)
            
        # Button to add a file
        self.add_button = tk.Button(self.analysis_frame, text='add analysis to stack', font=self.font_settings, width=int(self.maxwidth/4.3125), command=self.add2stack)
        self.add_button.grid(column=0, row=self.nanalysis+1, columnspan=1)
            
        # Button to remove a file
        self.remove_button = tk.Button(self.analysis_frame, text='remove last analysis from stack', font=self.font_settings, width=int(self.maxwidth/3), command=self.remove_last)
        self.remove_button.grid(column=1, row=self.nanalysis+1, sticky='news', columnspan=2)
        
        # Button to clear all files
        self.clear_button = tk.Button(self.analysis_frame, text='clear stack', font=self.font_settings, width=int(self.maxwidth/2.925), command=self.clear_all)
        self.clear_button.grid(column=3, row=self.nanalysis+1, columnspan=1)
        
        # Button load compute help
        self.analysis_help = tk.Button(self.analysis_frame, text='analysis help', font=self.font_settings, width=int(self.maxwidth/3.6), command=self.analysis_options)
        self.analysis_help.grid(column=4, row=self.nanalysis+1, columnspan=1)
            
        # Add padding to all frames in self.analysis_frame
        for widget in self.analysis_frame.winfo_children():
            widget.grid_configure(padx=int(self.xpadding/2), pady=int(self.ypadding/3))
            
        
        #-------------------------#
        # Initalize execute frame #
        #-------------------------#
        self.execute_frame = tk.LabelFrame(self.frame, borderwidth=0, highlightthickness=0)
        self.execute_frame.grid(column=0, row=4, padx=self.xpadding, pady=self.ypadding)
                
        # data2csv
        self.data2csv_btn = tk.Button(self.execute_frame, width=int(self.maxwidth/2.4), text='write loaded data to csv', command=self.write2csv, font=self.font_settings)
        self.data2csv_btn.grid(column=0, row=0)
        
        # save_mode
        self.save_mode_btn = tk.Button(self.execute_frame, width=int(self.maxwidth/2.4), text='save settings as mode', command=self.save_mode, font=self.font_settings)
        self.save_mode_btn.grid(column=1, row=0)
        
        # update plot
        self.update_plot_btn = tk.Button(self.execute_frame, width=int(self.maxwidth/2.4), text='update plot', command=self.analyze_and_plot, font=self.font_settings)
        self.update_plot_btn.grid(column=2, row=0)
        

        
        # Add padding to all frames in self.analysis_frame
        for widget in self.execute_frame.winfo_children():
            widget.grid_configure(padx=self.xpadding, pady=int(self.ypadding/3))
        
        #------------------------#
        # Run mainloop and close #
        #------------------------#
        self.root.protocol('WM_DELETE_WINDOW', self.closing)
        self.root.mainloop()
    
    #################################
    # Functions to call as commands #
    #################################
    # Function to save current GUI settings as a mode
    def save_mode(self):
        # Function to get directory
        def directory_path():
            path =filedialog.askdirectory(initialdir=self.pwd)
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
                    if len(self.methods) > 0 and self.methods.count('skip') < len(self.methods):
                        analysis = [] # [method, xlo, xhi, miscs, names]
                        for i, j, k, l, m  in zip(self.methods, self.xlos, self.xhis, self.miscs, self.names):
                            try: 
                                method = i.get(); xlo = j.get(); xhi = k.get(); misc = l.get(); name = m.get();
                                if method != 'skip':
                                    analysis.append([method, xlo, xhi, misc, name])
                            except: pass
                        f.write('analysis = [')
                        for n, (method, xlo, xhi, misc, name) in enumerate(analysis, 1):
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
                    f.write("{:>8}{}: {}\n".format('',"'analysis'", 'analysis'))
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
        
    # Function to toggle if replacing log file during loading mode
    def print_selection(self):
        if (self.replace.get()) == 1:
            self.log.out(f'Will replace logfile when loading mode (var = {self.replace.get()})')
        else:
            self.log.out(f'Will NOT replace logfile when loading mode (var = {self.replace.get()})')
            
    # Function to get filepath for modefile
    def modefile_path(self):
        ftypes = (('Python files (.py)', '*.py'), ('all files', '*.*'))
        startpath = io_functions.path_to_string(os.path.join(self.pwd, self.modespath))
        path = filedialog.askopenfilename(initialdir=startpath, title='Open mode file?', filetypes=ftypes)
        if path:
            self.modespath = os.path.dirname(os.path.abspath(path))
            path = io_functions.path_to_string(os.path.relpath(path))
            self.modefile.delete(0, tk.END); self.modefile.insert(0, path);
            
            module = misc_funcs.import_file(path)
            mode = module.mode
            self.load_mode(mode)
        else: self.log.GUI_error('ERROR could not load mode')
    
    # Function to load mode
    def load_mode(self, mode):       
        # Start updating settings
        if self.load_replace_logfile.get() == 'True':
            self.logfile.delete(0, tk.END)
            self.logfile.insert(0, mode['logfile'])
        
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
        
        # Start updating analysis
        self.clear_all()
        nloaded = len(self.methods)
        analysis = mode['analysis']
        for n, i in enumerate(analysis):
            method, xlo, xhi, misc, name = i
            if n >= nloaded: # add more boxes as needed
                self.add_overloaded_analysis()
            self.methods[n].current(self.supported_methods.index(method))
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
    
    # Function to analyze data and plot
    def analyze_and_plot(self):
        start_time = time.time()
        self.log.clear_all()
        self.log.out('\n\n--------------------------------------------------------------------------------------')
        self.log.out('                              Inputs and Data2load sections')
        self.log.out('--------------------------------------------------------------------------------------')
        # Read logfile and get data
        try:
            logfile = self.logfile.get()
            self.log.out('  logfile={}'.format(logfile))
            if os.path.isfile(logfile):
                keywords = self.keywords_entry.get().split(',')
                sections = self.sections_entry.get()
                self.log.out('  keywords={}'.format(self.keywords_entry.get()))
                self.log.out('  sections={}'.format(sections))
                log = read_log.file(logfile, keywords=keywords)
                data = log.get_data(sections, pflag=True) # {column-name:[lst of data]}
                self.columns = list(data.keys()) # update columns to push to xdata, ydata drop downs
            else: data = {}; self.log.GUI_error(f'ERROR lammps logfile {logfile} does not exist');
        except: data = {}; self.log.GUI_error('ERROR failed to read logfile and load data');
            
        # create plot
        if data:            
            # Get data
            xdata = self.xdata_dropdown.get()
            ydata = self.ydata_dropdown.get()
            self.log.out('  xdata={}'.format(xdata))
            self.log.out('  ydata={}'.format(ydata))
            
            # Perform compute on data
            xcompute = self.xcompute_entry.get()
            ycompute = self.ycompute_entry.get()
            if xcompute:
                x, message = misc_funcs.compute_thermo_data(data, xcompute)
                self.log.out('  xcompute={}'.format(xcompute))
                self.log.out('  {}'.format(message))
                xdata = 'xcompute'
            else: x = data[xdata] 
            if ycompute:
                y, message = misc_funcs.compute_thermo_data(data, ycompute)
                self.log.out('  ycompute={}'.format(ycompute))
                self.log.out('  {}'.format(message))
                ydata = 'ycompute'
            else: y = data[ydata] 
            
            # Log Labels
            self.log.out('  xlabel={}'.format(self.xlabel_entry.get()))
            self.log.out('  xlabel={}'.format(self.ylabel_entry.get()))
            
            # Intialize fig and try to close any currently open plots
            try: plt.close()
            except: pass
            
            # function to build data2plot dict. The following meaings:
            #    x = list/array of xdata to plot
            #    y = list/array of xdata to plot
            #    style = 'marker' or 'line' or 'both' or 'horizontal' or 'vertical'
            #    marker = marker style
            #    line = line style
            #    size = int (for 'point' or 'line' style)
            #    label = 'to put in legend'
            #    shiftable = Boolean, whether the data is shiftable or not (lin-reg for stress-v-strain)
            def plot_parms(x=[], y=[], style='point', marker='.', line='-', size=4, label='default', shiftable=False):
                plot_setup = {}
                plot_setup['x'] = x
                plot_setup['y'] = y
                plot_setup['size'] = size
                plot_setup['style'] = style
                plot_setup['marker'] = marker
                plot_setup['line'] = line
                plot_setup['label'] = label
                plot_setup['shiftable'] = shiftable
                return plot_setup
            
            # Function to format analysis for log file
            def format_analysis(method, xlo, xhi, misc, name):
                text = '\n\n--------------------------------------------------------------------------------------\n'
                txt = 'method={}; xlo={}; xhi={}; misc={}; name="{}"'.format(method, xlo, xhi, misc, name)
                chunks = len(text)
                nchunks = math.ceil(len(txt)/chunks)
                chunk_size = math.ceil(len(txt)/nchunks)
                tmp = [ txt[i:i+chunk_size] for i in range(0, chunks, chunk_size) ]
                for i in tmp:
                    text += '{}\n'.format(i)
                text += '--------------------------------------------------------------------------------------'
                return text
            
            # Get any anaylsis that users may want
            analysis = []; apply_moving_average = False; apply_butterworth_filter = False;
            LAMMPS_data_misc = ''; LAMMPS_data_xlo = ''; LAMMPS_data_xhi = '';
            for n, (i, j, k, l, m)  in enumerate(zip(self.methods, self.xlos, self.xhis, self.miscs, self.names)):
                try: 
                    method = i.get(); xlo = j.get(); xhi = k.get(); misc = l.get(); name = m.get();
                    if method == 'skip': continue
                    if method != '' and xlo != '':
                        try: xlo = float(xlo)
                        except: xlo = 'min-of-xdata'; self.log.GUI_error(f'ERROR xlo {xlo} is not a float, using minimum of xdata instead')
                    else: 
                        xlo = 'min-of-xdata'
                    if method != '' and xhi != '': 
                        try: xhi = float(xhi)
                        except: xhi = 'max-of-xdata'; self.log.GUI_error(f'ERROR xhi {xhi} is not a float, using maximum of xdata instead')
                    else: 
                        xhi = 'max-of-xdata'
                    if name == '' and method not in ['LAMMPS data (remove from plot)', 'LAMMPS data (apply moving average)', 'LAMMPS data (apply butterworth filter)']:
                        name = 'analysis-{}'.format(n)
                        self.log.warn(f'WARNING name was left empty, imposing {name}')
                    if method == 'LAMMPS data (apply moving average)':
                        apply_moving_average = True; LAMMPS_data_misc = misc;
                        LAMMPS_data_xlo = xlo; LAMMPS_data_xhi = xhi; continue
                    if method == 'LAMMPS data (apply butterworth filter)':
                        apply_butterworth_filter = True; LAMMPS_data_misc = misc;
                        LAMMPS_data_xlo = xlo; LAMMPS_data_xhi = xhi; continue
                    analysis.append([method, xlo, xhi, misc, name])
                except: pass

            
            # Save LAMMPS data to plot
            rm_lmp_data = False
            lmpdata = plot_parms(x=x, y=y, style='point', marker='.', line='-', size=4, label='LAMMPS data', shiftable=True)
            data2plot = [lmpdata]
            
            # Apply moving average or filter to LAMMPS data based on inputs
            if apply_moving_average and apply_butterworth_filter:
                self.log.GUI_error('ERROR can not use both "LAMMPS data (apply moving average)" and "LAMMPS data (apply butterworth filter)" at a time')
            else:
                try: LAMMPS_data_xlo = float(LAMMPS_data_xlo)
                except: LAMMPS_data_xlo = min(x)
                try: LAMMPS_data_xhi = float(LAMMPS_data_xhi)
                except: LAMMPS_data_xhi = max(x)
                
                if apply_moving_average:
                    LAMMPS_data_label = 'LAMMPS data w/moving average'
                    setting = self.get_misc_setting(LAMMPS_data_misc)
                    misc = LAMMPS_data_misc
                    if 'window' in LAMMPS_data_misc:
                        window = setting['window']
                    else: 
                        window = 100;
                        misc = ' default-window=100';
                    
                    self.log.out('  Implementing: LAMMPS data (apply moving average) with {} settings'.format(misc))
                    x, y = self.moving_average(x, y, LAMMPS_data_xlo, LAMMPS_data_xhi, window)
                    data2plot.append(plot_parms(x=x, y=y, style='point', marker='.', line='-', size=4, label=LAMMPS_data_label, shiftable=True))
                if apply_butterworth_filter:
                    LAMMPS_data_label = 'LAMMPS data w/butterworth filter'
                    setting = self.get_misc_setting(LAMMPS_data_misc)
                    misc = LAMMPS_data_misc
                    if 'order' in setting:
                        order = setting['order']
                    else: order=2; misc += ' default-order=2;';
                    
                    if 'wn' in setting:
                        wn = setting['wn']
                    else: wn=0.01; misc += ' default-wn=0.01;';
                    
                    self.log.out('  Implementing: LAMMPS data (apply butterworth filter) with {} settings'.format(misc))
                    x, y = self.butterworth_lowpass(x, y, LAMMPS_data_xlo, LAMMPS_data_xhi, wn, order)
                    data2plot.append(plot_parms(x=x, y=y, style='point', marker='.', line='-', size=4, label=LAMMPS_data_label, shiftable=True))
            
            # Start performing any analysis if present
            shift=False; shift_amount = 0; # Defaults
            if analysis:
                for method, xlo, xhi, misc, name in analysis:
                    if xlo == 'min-of-xdata':
                        xlo = min(x)
                        self.log.out('  xlo was not provided, using minimum of xdata {}'.format(xlo))
                    if xhi == 'max-of-xdata':
                        xhi = max(x)
                        self.log.out('  xhi was not provided, using maximum of xdata {}'.format(xhi))
                    
                    if method == 'LAMMPS data (remove from plot)':
                        rm_lmp_data = True
                    
                    if method == 'average':
                        reduced_x, reduced_y, average_y = self.average(x, y, xlo, xhi)
                        label = '{} (average = {:.6f})'.format(name, average_y)
                        self.log.out(format_analysis(method, xlo, xhi, misc, name))
                        self.log.out('  {}'.format(label))
                        avgdata = plot_parms(x=reduced_x, y=reduced_y, style='point', marker='.', line='-', size=4, label=label, shiftable=True)
                        data2plot.append(avgdata)
                        
                    if method == 'moving average':
                        setting = self.get_misc_setting(misc)
                        if 'window' in setting:
                            window = setting['window']
                        else: 
                            window = 100;
                            misc += ' default-window=100';
                        x_movavg, y_movavg = self.moving_average(x, y, xlo, xhi, window)
                        label = '{} ({})'.format(name, misc)
                        self.log.out(format_analysis(method, xlo, xhi, misc, name))
                        self.log.out('  {}'.format(label))
                        movavgdata = plot_parms(x=x_movavg, y=y_movavg, style='point', marker='.', line='-', size=4, label=label, shiftable=True)
                        data2plot.append(movavgdata)
                        
                    if method == 'butterworth (low pass)':
                        if 'order' in misc or 'cutoff':
                            setting = self.get_misc_setting(misc)
                            if 'order' in setting:
                                order = setting['order']
                            else: order=2; misc += ' default-order=2;';
                            
                            if 'wn' in setting:
                                wn = setting['wn']
                            else: wn=0.01; misc += ' default-wn=0.01;';

                        label = '{} wn = {}; cutoff = {}'.format(name, order, wn)
                        self.log.out('  {}'.format(label))
                        xfilter, yfilter = self.butterworth_lowpass(x, y, xlo, xhi, wn, order)
                        self.log.out(format_analysis(method, xlo, xhi, misc, name))
                        bwfdata = plot_parms(x=xfilter, y=yfilter, style='point', marker='.', line='-', size=4, label=label, shiftable=True)
                        data2plot.append(bwfdata)
                    
                    if method == 'linear regression':
                        if 'shift' in misc or 'extend':
                            setting = self.get_misc_setting(misc)
                            if 'shift' in setting:
                                shift = setting['shift']
                            else: shift=False; misc = 'default-shift=False';
                        else: 
                            shift=False
                            misc += ' default-shift=False'
                        b0, b1, xreg, yreg = self.linear_regression(x, y, xlo, xhi)
                        self.log.out(format_analysis(method, xlo, xhi, misc, name))
                        if shift: 
                            if b0 < 0:
                                shift_amount = abs(b0)
                            else: shift_amount = -abs(b0)
                            self.log.out(f'  Shifting all Y-data by {-shift_amount}')
                            b0 = 0
                        label = '{} y = {:.6f}x + {:.6f} (shifted by {:.2f})'.format(name, b1, b0, -shift_amount)
                        self.log.out('  {}'.format(label))
                        if 'extend' in setting:
                            xreg = [xlo, xhi]
                            extend = setting['extend']
                            if extend > 0:
                                xreg.insert(2, xreg[1]+extend)
                            else: xreg.insert(0, xreg[0]+extend)
                            yreg = [i * b1 + b0 for i in xreg]
                            regdata = plot_parms(x=xreg, y=yreg, style='both', marker='o', line='-', size=3, label=label, shiftable=True)
                        else:
                            regdata = plot_parms(x=xreg, y=yreg, style='line', marker='.', line='-', size=3, label=label, shiftable=True)
                        data2plot.append(regdata)
                        
                    if method in ['minimum', 'maximum']:
                        setting = self.get_misc_setting(misc)
                        if 'window' in setting:
                            window = setting['window']
                        else: 
                            window = 100;
                            misc += ' default-window=100';
                        x_movavg, y_movavg = self.moving_average(x, y, xlo, xhi, window)
                        label = '{} ({})'.format(name, misc)
                        movavgdata = plot_parms(x=x_movavg, y=y_movavg, style='point', marker='.', line='-', size=4, label=label, shiftable=True)
                        data2plot.append(movavgdata)
                        self.log.out(format_analysis(method, xlo, xhi, misc, name))
                        self.log.out('  {}'.format(label))
                        
                        x_movavg = list(x_movavg); y_movavg = list(y_movavg)
                        if method == 'minimum': yy = min(y_movavg)
                        if method == 'maximum': yy = max(y_movavg)
                        xx =  x_movavg[y_movavg.index(yy)]
                        label = '{} (x={}, y={})'.format(name, xx, yy)
                        minmaxdata = plot_parms(x=xx, y=yy, style='point', marker='o', line='-', size=8, label=label, shiftable=True)
                        data2plot.append(minmaxdata)
                        self.log.out('  {}'.format(label))
                    
                    if method == 'cursor':
                        setting = self.get_misc_setting(misc)
                        if 'x' in setting or 'y' in setting:
                            if 'x' in setting:
                                xvalue = setting['x']
                            else: xvalue = None
                            
                            if 'y' in setting:
                                yvalue = setting['y']
                            else: yvalue = None
                            
                            # Set style and marker based on if x or y or x & y are set
                            plot_and_log = True
                            if xvalue is not None and yvalue is not None:
                                cursor_style = 'point'
                                cursor_marker = '+'
                            elif xvalue is not None:
                                cursor_style = 'vertical'
                                cursor_marker = ':'
                            elif yvalue is not None:
                                cursor_style = 'horizontal'
                                cursor_marker = ':'
                            else: plot_and_log = False
                            
                            if plot_and_log:
                                label = '{} ({})'.format(name, misc)
                                self.log.out(format_analysis(method, xlo, xhi, misc, name))
                                self.log.out('  {}'.format(label))
                                cursordata = plot_parms(x=xvalue, y=yvalue, style=cursor_style, marker=cursor_marker, line='-', size=4, label=label, shiftable=False)
                                data2plot.append(cursordata)
                    
                    if method == 'hyperbola':
                        minimum_convergence=None
                        initial_guess=False
                        setting = self.get_misc_setting(misc)
                        if 'p' in setting:
                            minimum_convergence = setting['p']
                        if 'initial_guess' in setting:
                            initial_guess = setting['initial_guess']
                        xout, yout, params, center, slopes, transition = self.fit_hyperbola(x, y, xlo, xhi, minimum_convergence=minimum_convergence, initial_guess=initial_guess, maxiter=10**6)
                        self.log.out(format_analysis(method, xlo, xhi, misc, name))
                        
                        label = '{} slopes (lower-slope={};\nupper-slope={})'.format(name, slopes[0], slopes[1])
                        self.log.out('  {} - slopes (lower-slope={}; upper-slope={})'.format(name, slopes[0], slopes[1]))
                        hyperboladata = plot_parms(x=xout, y=yout, style='point', marker='.', line='-', size=2, label=label, shiftable=False)
                        data2plot.append(hyperboladata)
                        
                        label = '{} center (x={:.4f}; y={:.4f})'.format(name, center[0], center[1])
                        self.log.out('  {}'.format(label))
                        centerdata = plot_parms(x=center[0], y=center[1], style='point', marker='.', line='-', size=10, label=label, shiftable=False)
                        data2plot.append(centerdata)
                        
                        if minimum_convergence is not None:
                            label = '{} transition region (P={}; xlo={:.4f}; xhi={:.4f})'.format(name, minimum_convergence, transition[0], transition[1])
                            self.log.out('  {}'.format(label))
                            miny = [min(yout) for _ in transition]
                            transitiondata = plot_parms(x=transition, y=miny, style='both', marker='|', line='-', size=6, label=label, shiftable=False)
                            data2plot.append(transitiondata)
                            
                    if method == 'piecewise-regression':
                        setting = self.get_misc_setting(misc)
                        if 'n' in setting:
                            n = setting['n']
                        else: n = 1; misc += ' default-n=1';
                        if 'shift' in setting:
                            shift = setting['shift']
                        else: shift = False; misc += ' default-shift=False';
                            
                        xout, yout, xbreaks, ybreaks, slopes = self.piecewise_regression(x, y, xlo, xhi, n)
                        self.log.out(format_analysis(method, xlo, xhi, misc, name))
                        
                        shift_amount = 0
                        if shift:
                            b0 = ybreaks[0]
                            if b0 < 0:
                                shift_amount = abs(b0)
                            else: shift_amount = -abs(b0)
                            self.log.out(f'  Shifting all Y-data by {-shift_amount}')
                        
                        label = '{} (n-breakpoints={})'.format(name, n)
                        self.log.out(f'  {label}')
                        piecewisedata = plot_parms(x=xout, y=yout, style='line', marker='.', line='-', size=3, label=label, shiftable=True)
                        data2plot.append(piecewisedata)
                        
                        label = '{} slopes and points: \n'.format(name)
                        for i, j in slopes:
                            slope = slopes[(i, j)]
                            label += '   slope between (p{}, p{}) = {}\n'.format(i, j, slope)
                            label += '    p{} = ({}, {})\n'.format(i, xbreaks[i], ybreaks[i]+shift_amount)
                            label += '    p{} = ({}, {})\n'.format(j, xbreaks[j], ybreaks[j]+shift_amount)
                        self.log.out(f'  {label}')
                        breakdata = plot_parms(x=xbreaks, y=ybreaks, style='point', marker='.', line='-', size=10, label=label, shiftable=True)
                        data2plot.append(breakdata)
                        
                    if method == 'spline-integration':
                        setting = self.get_misc_setting(misc)
                        if 'window' in setting:
                            window = setting['window']
                        else: 
                            window = 100;
                            misc += '  default-window=100';
                        if 'shift' in setting:
                            shift = setting['shift']
                        else: shift=False; misc += '  default-shift=False';
                            
                        self.log.out(format_analysis(method, xlo, xhi, misc, name))
                        xout, yout, area, shift_amount = self.spline_integration(x, y, xlo, xhi, window, shift)
                        label = '{} (area: {})'.format(name, area)
                        self.log.out(f'  {label}')
                        inetgrateddata = plot_parms(x=xout, y=yout, style='point', marker='.', line='-', size=6, label=label, shiftable=False)
                        data2plot.append(inetgrateddata)
                            
                        
            # Generate plot
            fig, ax = plt.subplots()
            colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
            color_index = 0
            for data in data2plot:
                xx = data['x']
                yy = data['y']
                label = data['label']
                style = data['style']
                shiftable = data['shiftable']
                marker = data['marker']
                line = data['line']
                size = data['size']
                color = colors[color_index]
                
                # Skip lammps data if rm_lmp_data and label is 'LAMMPS data'
                if rm_lmp_data and label == 'LAMMPS data': continue
                
                # If shift, shift all data
                if shift and shiftable:
                    try: yy = [i+shift_amount for i in yy]
                    except: 
                        try: yy += shift_amount
                        except: pass
                
                # Plot data with or without line width size
                if style == 'line':
                    plt.plot(xx, yy, line, color=color, lw=size, label=label)
                elif style == 'point':
                    plt.plot(xx, yy, marker, color=color, ms=size, label=label)
                elif style == 'both':
                    plt.plot(xx, yy, marker, color=color, ls='-', ms=size, label=label)
                elif style == 'vertical':
                    plt.axvline(xx, color=color, ls=marker, lw=size, label=label)
                elif style == 'horizontal':
                    plt.axhline(yy, color=color, ls=marker, lw=size, label=label)
                else: raise Exception(f'ERROR {style} not supported')
                
                # increment color index and rest to zero if to large
                color_index += 1
                if color_index + 1 > len(colors):
                    color_index = 0
            
            # Set labels
            if self.xlabel_entry.get():
                ax.set_xlabel(self.xlabel_entry.get())#, fontsize=int(2*self.font_size))
            else: ax.set_xlabel(xdata)
            if self.ylabel_entry.get():
                ax.set_ylabel(self.ylabel_entry.get())#, fontsize=int(2*self.font_size))
            else: ax.set_ylabel(ydata)

            
            # Set legend and size
            ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fancybox=True, ncol=1)#, fontsize=int(2*self.font_size))
            fig.set_size_inches(12, 4, forward=True)
            fig.tight_layout()
            plt.show()
            
            # Save current image
            if self.settings['save-fig']:
                figname = '{}_X={}_Y={}.jpeg'.format(self.logfile.get(), xdata, ydata)
                fig.savefig(figname, dpi=self.settings['image-dpi'])
        else: self.log.GUI_error('ERROR no data loaded, unable to plot none existant data')
        
        # Script run time
        execution_time = (time.time() - start_time)
        self.log.out('Execution time in seconds: ' + str(execution_time))
        
        # write log
        logname = '{}_X={}_Y={}.log.lunar'.format(self.logfile.get(), xdata, ydata)
        self.log.write_logged(logname)
        #self.popup(self.log.logged, title='Outputs', width=150)
        return
    
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
    
    ###############################################
    # Function to get misc dict of key/value pair #
    ###############################################
    def get_misc_setting(self, misc):
        setting = {} # {keyword:float or int or Boolean}
        tmp1 = misc.split(';')
        for tmp2 in tmp1:
            tmp3 = tmp2.split('=')
            if len(tmp3) >= 2:
                i = tmp3[0].strip()
                try: j = eval(tmp3[1])
                except: j = str(tmp3[1])
                setting[i] = j
        return setting
    
    #############################
    # average analysis function #
    #############################
    def average(self, x, y, xlo, xhi):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            average_y = misc_funcs.avg(reduced_y)
        else:
            average_y = 0
            self.log.GUI_error(f'ERROR (average) no LAMMPS data in xrange {xlo} - {xhi}')
        return reduced_x, reduced_y, average_y
    
    #########################################################
    # The method below implements a moving average function #
    #########################################################
    def moving_average(self, x, y, xlo, xhi, window):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            x_movavg = misc_funcs.moving_average(reduced_x, window)
            y_movavg = misc_funcs.moving_average(reduced_y, window)
        else:
            x_movavg = [0, 1]; y_movavg = [0, 1]
            self.log.GUI_error(f'ERROR (moving average) no LAMMPS data in xrange {xlo} - {xhi}')
        return list(x_movavg), list(y_movavg)
    
    ############################################################
    # The method below implements a lowpass butterworth filter #
    ############################################################
    def butterworth_lowpass(self, x, y, xlo, xhi, cutoff, order):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            yfilter = misc_funcs.butter_lowpass_filter(reduced_x, reduced_y, cutoff, order)
        else:
            yfilter = [i for i in reduced_y]
            self.log.GUI_error(f'ERROR (butterworth low pass) no LAMMPS data in xrange {xlo} - {xhi}')
        return list(reduced_x), list(yfilter)
    
    
    #########################################################
    # The method below implements a linear regression model #
    #########################################################
    def linear_regression(self, x, y, xlo, xhi):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            linear_regression = misc_funcs.linear_regression(reduced_x, reduced_y)
            b1 = linear_regression.b1; b0 = linear_regression.b0; xreg = [xlo, xhi];
            yreg = [i * b1 + b0 for i in xreg]
        else:
            xreg = [0, 1]; yreg = [0, 1]
            self.log.GUI_error(f'ERROR (linear regression) no LAMMPS data in xrange {xlo} - {xhi}')
        return b0, b1, xreg, yreg
    
    ####################################################
    # The method below implements a spline-integration #
    ####################################################
    def spline_integration(self, x, y, xlo, xhi, window, shift):
        from scipy.interpolate import InterpolatedUnivariateSpline
        x_movavg, y_movavg = self.moving_average(x, y, xlo, xhi, window)
        if x_movavg.size != 0 and y_movavg.size != 0:
            xs = np.array(x_movavg)
            ys = np.array(y_movavg)
            if shift:
                if ys[0] < 0:
                    shift_amount = abs(ys[0])
                else: shift_amount = -abs(ys[0])   
                ys = ys + shift_amount               
            else: shift_amount = 0
            
            spl = InterpolatedUnivariateSpline(xs, ys, k=1)  # k=1 gives linear interpolation
            area = spl.integral(xlo, xhi)
            xout = list(xs)
            yout = list(spl(xs))
        else:
            xout = [0, 1]; yout = [0, 1]; area = 0;
            self.log.GUI_error(f'ERROR (spline-integration) no LAMMPS data in xrange {xlo} - {xhi}')
        return xout, yout, area, shift_amount
    
    ####################################################################
    # The method below implements a "piecewise regression" on the data #
    ####################################################################
    def piecewise_regression(self, x, y, xlo, xhi, n):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:
            xout, yout, xbreaks, ybreaks, slopes = misc_funcs.piecewise_regression(reduced_x, reduced_y, xlo, xhi, n)
        else:
            xout = [0, 1]; yout = [0, 1]; slopes = {(0,1):1}
            xbreaks = [0, 1]; ybreaks = [0, 1];
            self.log.GUI_error(f'ERROR (peicewise-regression) no LAMMPS data in xrange {xlo} - {xhi}')
        return xout, yout, xbreaks, ybreaks, slopes
    
    ######################################################################################
    # The method below implements the Tg/CTE hyperbola fit from the following paper:     #
    # Uncertainty quantification in molecular dynamics studies of the glass transition   #
    # temperature - Paul N. Patrone, Andrew Dienstfrey, ... - Polymer Volume 87 - 2016   #
    ######################################################################################
    def fit_hyperbola(self, x, y, xlo, xhi, minimum_convergence=None, initial_guess=False, maxiter=10**4):
        reduced_x, reduced_y = misc_funcs.reduce_data(x, y, xlo, xhi)
        if reduced_x and reduced_y:        
            xout, yout, params, center, slopes, transition = misc_funcs.fit_hyperbola(x, y, xlo, xhi, minimum_convergence, initial_guess, maxiter)
        else:
            xout = [0, 1]; yout = [0, 1]; slopes = [0, 1]
            center = [0, 0]; params = [0, 0, 0, 0, 0];
            transition = [];
            self.log.GUI_error(f'ERROR no (hyperbola) LAMMPS data in xrange {xlo} - {xhi}')
        return xout, yout, params, center, slopes, transition
    
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
        
        self.xlo = tk.Entry(self.analysis_frame, width=int(self.maxwidth/5), font=self.font_settings)
        self.xlo.grid(column=1, row=self.nanalysis)
        #self.xlo.insert(0, 0)
        self.xlos.append(self.xlo)
        
        self.xhi = tk.Entry(self.analysis_frame, width=int(self.maxwidth/5), font=self.font_settings)
        self.xhi.grid(column=2, row=self.nanalysis)
        #self.xhi.insert(0, 1)
        self.xhis.append(self.xhi)
        
        self.misc = tk.Entry(self.analysis_frame, width=int(self.maxwidth/2), font=self.font_settings)
        self.misc.grid(column=3, row=self.nanalysis)
        self.misc.insert(0, '')
        self.miscs.append(self.misc)
        
        self.name = tk.Entry(self.analysis_frame, width=int(self.maxwidth/3), font=self.font_settings)
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
            widget.grid_configure(padx=int(xpadding/2), pady=int(ypadding/3))
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
        path = filedialog.askopenfilename(initialdir=self.filepath, title='Open logfile?', filetypes=ftypes)
        if path:
            self.filepath = os.path.dirname(os.path.abspath(path))
            path = os.path.relpath(path)
            self.logfile.delete(0, tk.END); self.logfile.insert(0, path);
            
            # Update self.columns
            if os.path.isfile(path):
                try: keywords = self.keywords.get().split(',')
                except: keywords = ['Step', 'Temp', 'Density', 'Press', 'TotEng', 'KinEng', 'PotEng', 'Lx', 'Ly', 'Lz']
                try: 
                    sections = self.sections_entry.get()
                    if sections == '': sections = 'all'
                except: sections = 'all'
                log = read_log.file(path, keywords=keywords)
                data = log.get_data(sections, pflag=True) # {column-name:[lst of data]}
                self.columns = list(data.keys())
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
                data2write = log.get_data(sections, pflag=True) # {column-name:[lst of data]}
                self.columns = list(data2write.keys()) # update columns to push to xdata, ydata drop downs
            else: data2write = {}; print(f'lammps logfile {logfile} does not exist');
        except: data2write = {}; print('ERROR failed to read logfile and load data')
        
        # writedata to csv
        if data2write:
            csvname = '{}.csv'.format(self.logfile.get())
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