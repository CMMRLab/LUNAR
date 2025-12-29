# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
February 20, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.GUI_scale_settings as GUI_scale_settings
import src.atom_typing.log_processor.main as main
import src.io_functions as io_functions
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
            if isinstance(settings['mode'], str):
                mode = main.import_file(settings['mode']).mode
            else: mode = settings['mode']
        except:
            # dependant and independant variable strings
            independent_variable = "numeric"
            dependent_variables = """
            rings[6]['%Ring'],    rings[6]['C']['%natoms'],    hybridizations['Sp1-C']['%natoms']
            hybridizations['Sp2-C']['%natoms'],    hybridizations['Sp3-C']['%natoms'],    hybridizations['all-H']['%natoms']
            """
            
            # logger string
            logger = """
            filename,    numeric,    wildcards,    rings[6]['%Ring'],    rings[6]['C']['%natoms'],    hybridizations['Sp1-C']['%natoms']
            hybridizations['Sp2-C']['%natoms'],    hybridizations['Sp3-C']['%natoms'],    hybridizations['all-H']['%natoms']
            """
            # loadable mode
            mode = {'independent_variable': independent_variable,
                    'dependent_variables': dependent_variables,
                    'parent_directory': 'logfile',
                    'sorting_direction': 'ascending',
                    'sorting_method': 'sort',
                    'logfile': 'EXAMPLES/array_processing/atom_typing_logfile_processor/poly_tracking_replicate_1_time_*ps_typed.log.lunar',
                    'newfile': 'atom_typing_logfile_processor',
                    'logger': logger
                    }
            self.log.GUI_error(f'ERROR loading mode file {settings["mode"]}. Internally deriving default settings. Likely launching from outside of LUNAR directory.')
        self.mode = mode

        
        # Initialize window
        self.root = tk.Tk()
        self.root.title('LUNAR/atom_typing_log_processor.py GUI v1.0')
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
        self.inputs_frame.grid(row=0, column=0, sticky='news', padx=self.xpadding, pady=self.ypadding)
        
        # logfile selection button
        self.logfile = tk.Entry(self.inputs_frame, width=int(1.495*self.maxwidth), font=self.font_settings)
        self.logfile.insert(0, self.mode['logfile'])
        self.logfile.grid(column=1, row=0, columnspan=4)
        self.logfile_button = tk.Button(self.inputs_frame, text='logfile', font=self.font_settings, command=self.logfile_path)
        self.logfile_button.grid(column=0, row=0)
                
        # parent_directory selection button
        self.parent_directory = tk.Entry(self.inputs_frame, width=int(1.495*self.maxwidth), font=self.font_settings)
        self.parent_directory.insert(0, self.mode['parent_directory'])
        self.parent_directory.grid(column=1, row=1, columnspan=4)
        self.parent_directory_button = tk.Button(self.inputs_frame, text='parent_directory', font=self.font_settings, command=self.directory_path)
        self.parent_directory_button.grid(column=0, row=1)
        
        # modes
        self.modefile = tk.Entry(self.inputs_frame, width=int(0.95*self.maxwidth), font=self.font_settings)
        if isinstance(settings['mode'], str): 
            modefile = settings['mode']
        else: modefile = 'input_type=dict'
        self.modefile.insert(0, modefile)
        self.modefile.grid(column=1, row=2)
        self.modefile_button = tk.Button(self.inputs_frame, text='mode file', font=self.font_settings, command=self.modefile_path)
        self.modefile_button.grid(column=0, row=2)        
        
        # load_replace_logfile drop down menu
        styles = [True, False]
        self.load_replace_logfile = ttk.Combobox(self.inputs_frame, values=styles, width=int(self.maxwidth/9.25), font=self.font_settings)
        self.load_replace_logfile.current(styles.index(settings['replace_logfile_when_loading_mode']))
        self.load_replace_logfile.grid(column=4, row=2)
        self.load_replace_logfile_label = tk.Label(self.inputs_frame, text='Replace logfile when loading mode', font=self.font_settings)
        self.load_replace_logfile_label.grid(column=3, row=2)
        
        # newfile entry
        self.newfile = tk.Entry(self.inputs_frame, width=int(1.495*self.maxwidth), font=self.font_settings)
        self.newfile.insert(0, self.mode['newfile'])
        self.newfile.grid(column=1, row=3, columnspan=4)
        self.newfile_label = tk.Label(self.inputs_frame, text='newfile', font=self.font_settings)
        self.newfile_label.grid(column=0, row=3)

        
        # Add padding to all frames in self.inputs_frame
        for widget in self.inputs_frame.winfo_children():
            widget.grid_configure(padx=int(self.xpadding/2), pady=int(self.ypadding/3))
            
        #--------------#
        # Logger frame #
        #--------------#
        # Initalize load data frame
        self.logger_frame = tk.LabelFrame(self.frame, text='logger', font=self.font_settings)
        self.logger_frame.grid(row=1, column=0, sticky='news', padx=self.xpadding, pady=self.ypadding)
        
        # logger entry
        self.logger_text = tk.Text(self.logger_frame, width=int(1.7*self.maxwidth), height=10, font=self.font_settings)
        self.logger_text.insert(tk.END, self.mode['logger'].strip())
        self.logger_text.grid(column=0, row=1, columnspan=1)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.logger_frame.winfo_children():
            widget.grid_configure(padx=int(self.xpadding/2), pady=int(self.ypadding/3))
            
        #----------------#
        # Analysis frame #
        #----------------#
        # Initalize analysis frame
        self.analysis_frame = tk.LabelFrame(self.frame, text='Analysis options', font=self.font_settings)
        self.analysis_frame.grid(row=3, column=0, sticky='news', padx=self.xpadding, pady=self.ypadding)
        
        # independent_variable entry
        self.independent_variable = tk.Entry(self.analysis_frame, width=int(0.75*self.maxwidth), font=self.font_settings)
        self.independent_variable.insert(0, self.mode['independent_variable'])
        self.independent_variable.grid(column=1, row=0, columnspan=1)
        self.independent_variable_label = tk.Label(self.analysis_frame, text='independent_variable', font=self.font_settings)
        self.independent_variable_label.grid(column=0, row=0)
        
        # dependent_variables text
        self.dependent_variables = tk.Text(self.analysis_frame, width=int(1.5*self.maxwidth), height=6, font=self.font_settings)
        self.dependent_variables.insert(tk.END, self.mode['dependent_variables'].strip())
        self.dependent_variables.grid(column=1, row=2, columnspan=5)
        self.dependent_variables_label = tk.Label(self.analysis_frame, text='dependent_variables', font=self.font_settings)
        self.dependent_variables_label.grid(column=0, row=2)
        
        # sorting method drop down
        self.supported_sorting_methods = ['natsort', 'sort', 'numericsort', 'wildcards[i]']
        try:
            current = self.supported_sorting_methods.index(self.mode['sorting_method'])
        except:
            self.supported_sorting_methods.append(self.mode['sorting_method'])
            current = self.supported_sorting_methods.index(self.mode['sorting_method'])
        self.sorting_method = ttk.Combobox(self.analysis_frame, values=self.supported_sorting_methods, width=int(self.maxwidth/7.5), font=self.font_settings)
        self.sorting_method.current(current)
        self.sorting_method.grid(column=3, row=0)
        self.sorting_method_label = tk.Label(self.analysis_frame, text='sorting_method', font=self.font_settings)
        self.sorting_method_label.grid(column=2, row=0)
        
        # sorting method drop down
        self.supported_sorting_directions = ['ascending', 'descending']
        self.sorting_direction = ttk.Combobox(self.analysis_frame, values=self.supported_sorting_directions , width=int(self.maxwidth/7.5), font=self.font_settings)
        self.sorting_direction.current(self.supported_sorting_directions.index(self.mode['sorting_direction']))
        self.sorting_direction.grid(column=5, row=0)
        self.sorting_direction_label = tk.Label(self.analysis_frame, text='sorting_direction', font=self.font_settings)
        self.sorting_direction_label.grid(column=4, row=0)
        
        # Add padding to all frames in self.analysis_frame
        for widget in self.analysis_frame.winfo_children():
            widget.grid_configure(padx=int(self.xpadding/4), pady=int(self.ypadding/3))
            
            
        #-------------------------#
        # Initalize execute frame #
        #-------------------------#
        self.execute_frame = tk.LabelFrame(self.frame, borderwidth=0, highlightthickness=0)
        self.execute_frame.grid(column=0, row=4, padx=self.xpadding, pady=self.ypadding)
        
        # run button
        self.run_btn = tk.Button(self.execute_frame, width=int(1.45*self.maxwidth), text='Run LUNAR/atom_typing_log_processor', command=self.run_main, font=self.font_settings)
        self.run_btn.grid(column=0, row=0, columnspan=2)
        
        # save_mode
        self.save_mode_btn = tk.Button(self.execute_frame, width=int(0.65*self.maxwidth), text='save settings as mode', command=self.save_mode, font=self.font_settings)
        self.save_mode_btn.grid(column=0, row=1, sticky='news')
        
        # Quick help
        self.quick_help = tk.Button(self.execute_frame, width=int(0.65*self.maxwidth), text='Quick help', font=self.font_settings, command=self.quickhelp)
        self.quick_help.grid(column=1, row=1, sticky='news')
        
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
    # Function to get filepath for logfile
    def logfile_path(self):
        ftypes = (('LUNAR log files (.lunar, .log, .txt)', '*.lunar *.log *.txt'), ('all files', '*.*'))
        path = filedialog.askopenfilename(title='Open logfile?', filetypes=ftypes)
        if path:
            path = os.path.relpath(path)
            self.logfile.delete(0, tk.END); self.logfile.insert(0, path);
        return
    
    # Quick help button
    def quickhelp(self):
        try: # Try to get text from GUI_help_page.txt file
            txt = os.path.join(self.pwd, 'src/GUI_quick_help_pages/atom_typing_log_processor.txt')
            logged = []
            with open(txt, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    # Strip comment's and split by whitespace
                    line = line.split('#')[0]
                    line = line.rstrip()
                    logged.append(line)
        except: # except something failed
            logged.append('FAILED to read LUNAR/src/GUI_quick_help_pages/atom_typing_log_processor.txt document.')
            logged.append('Most likely cause is the atom_typing.txt file was renamed in the')
            logged.append('LUNAR/src/GUI_quick_help_pages/ directory or directory names were changed.')
        self.popup(logged, title='Quick help')
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
    
    # Function to get directory
    def directory_path(self):
        path =filedialog.askdirectory()
        if path:
            path = os.path.relpath(path)
            self.parent_directory.delete(0, tk.END); self.parent_directory.insert(0, path);
        return
    

        
    # Closing command    
    def closing(self):
        print('Terminating log_analysis GUI'); self.root.destroy();
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
        modes_dir = self.modespath
        if os.path.isdir(modes_dir):
            path = filedialog.asksaveasfilename(defaultextension=".py",
                                                filetypes=[("Pyton files", "*.py")],
                                                initialdir=modes_dir)
        else:
            path = filedialog.asksaveasfilename(defaultextension=".py",
                                                filetypes=[("Pyton files", "*.py")])
            
        if path:
            # Generate mode dictionary from GUI options
            mode = self.generate_mode()
            with open(path, 'w') as f:
                f.write('# -*- coding: utf-8 -*-\n')
                f.write('"""\n')
                f.write('Created by atom_typing_log_processor.py to store a mode dictionary\n')
                f.write('called "mode", which will then allow all settings to be\n')
                f.write('loaded by clicking on this file within log_anaylsis.py.\n\n')
                #f.write('{}'.format(about.get('1.0', tk.END)))
                f.write('"""\n')
                
                f.write('\n# dependant and independant variable strings\n')
                f.write('independent_variable = "{}"\n'.format(mode['independent_variable']))
                f.write('dependent_variables = {}\n'.format('"""'))
                f.write('{}\n'.format(mode['dependent_variables']))
                f.write('{}\n'.format('"""'))
                
                f.write('\n# logger string\n')
                f.write('logger = {}\n'.format('"""'))
                f.write('{}\n'.format(mode['logger']))
                f.write('{}\n'.format('"""'))
                
                f.write('\n# loadable mode\n')
                f.write("mode = {")
                f.write("{}{}: {},\n".format('',"'independent_variable'", 'independent_variable'))
                f.write("{:>8}{}: {},\n".format('',"'dependent_variables'", 'dependent_variables'))
                f.write("{:>8}{}: '{}',\n".format('',"'parent_directory'", io_functions.path_to_string(mode['parent_directory'])))
                f.write("{:>8}{}: '{}',\n".format('',"'sorting_direction'", mode['sorting_direction']))
                f.write("{:>8}{}: '{}',\n".format('',"'sorting_method'", mode['sorting_method']))
                f.write("{:>8}{}: '{}',\n".format('', "'logfile'", io_functions.path_to_string(mode['logfile'])))
                f.write("{:>8}{}: '{}',\n".format('',"'newfile'", mode['newfile']))
                f.write("{:>8}{}: {}\n".format('',"'logger'", 'logger'))
                f.write("{:>8}{}\n\n".format('', '}'))
        
        self.log.out('Saving current GUI settings as a mode')
        return
    
    # Function to load mode
    def load_mode(self, mode):       
        # Start updating settings
        if self.load_replace_logfile.get() == 'True':
            self.logfile.delete(0, tk.END)
            self.logfile.insert(0, mode['logfile'])
            
        self.parent_directory.delete(0, tk.END)
        self.parent_directory.insert(0, mode['parent_directory'])
        
        self.newfile.delete(0, tk.END)
        self.newfile.insert(0, mode['newfile'])
    
        self.logger_text.delete('1.0', tk.END)
        self.logger_text.insert(tk.END, mode['logger'].strip())
        
        self.independent_variable.delete(0, tk.END)
        self.independent_variable.insert(0, mode['independent_variable'])
        
        self.dependent_variables.delete('1.0', tk.END)
        self.dependent_variables.insert(tk.END, mode['dependent_variables'].strip())
        
        try: method = mode['sorting_method']
        except: method = 'sorted'
        try:
            current = self.supported_sorting_methods.index(method)
        except:
            self.supported_sorting_methods.append(method)
            current = self.supported_sorting_methods.index(method)
        self.sorting_method.current(current)
        
        try: direction = mode['sorting_direction']
        except: direction = 'ascending'
        self.sorting_direction.current(self.supported_sorting_directions.index(direction))
        return
    
    # Method to generate a mode dictionary
    def generate_mode(self):
        mode = {}
        mode['independent_variable'] = self.independent_variable.get()
        mode['dependent_variables'] = self.dependent_variables.get(1.0, 'end-1c').strip()
        mode['parent_directory'] = self.parent_directory.get()
        mode['sorting_direction'] = self.sorting_direction.get()
        mode['sorting_method'] = self.sorting_method.get()
        mode['logfile'] = self.logfile.get()
        mode['newfile'] = self.newfile.get()
        mode['logger'] = self.logger_text.get(1.0, 'end-1c').strip()
        return mode
    
    
    ########################################################################
    # Function to analyze data and plot (this is the "heart" of this code) #
    ########################################################################
    def run_main(self):
        # Generate mode dictionary from GUI options
        mode = self.generate_mode()
        try:
            main.main(mode, log=self.log)
        except Exception:
            self.log.GUI_error(traceback.format_exc())
        self.popup(self.log.logged, title='Outputs')
        self.log.clear_all()
        return