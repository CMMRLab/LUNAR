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
import src.GUI_scale_settings as GUI_scale_settings
from src.auto_morse_bond.main import main
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


#############################
# LUNAR/auto_morse_bond GUI #
#############################
class auto_morse_bond_GUI:
    def __init__(self, topofile, morsefile, parent_directory, newfile, mass_map, min_bond_length, coeffs2skip,
                 radius_specs, alpha_specs, alpha_scale, files2write, atom_style, zero_effected_xterms,
                 bondbreak_scale, ff_class, include_type_labels, class2xe_update, include_rcut, GUI_zoom):
        
        # Pass certain inputs as attribute (These will only be able to be adjusted from the python script)
        self.mass_map = mass_map
        
        # Find present working directory
        self.pwd = os.getcwd()
        self.filename = os.path.join(self.pwd, 'auto_morse_bond_update.py')

        # Initialize window
        self.root = tk.Tk()
        self.root.title('LUNAR/auto_morse_bond_update.py GUI v1.0')
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
        xpadding = int(math.ceil(GUI_SF*self.xpadding))
        ypadding = int(math.ceil(GUI_SF*self.ypadding))
        maxwidth = int(math.ceil(GUI_SF*self.maxwidth))
        font_settings = (self.font_type, font_size)
        
        #--------------#
        # Inputs frame #
        #--------------#
        # Initalize  inputs frame
        self.inputs_frame = tk.LabelFrame(self.frame1, text='Inputs', font=font_settings)
        self.inputs_frame.grid(row=0, column=0, columnspan=2, padx=xpadding, pady=ypadding)
        
        # topofile selection button
        self.topofile = tk.Entry(self.inputs_frame, width=int(1.15*maxwidth), font=font_settings)
        self.topofile.insert(0, topofile)
        self.topofile.grid(column=1, row=0)
        self.topofile_button = tk.Button(self.inputs_frame, text='topofile', font=font_settings, command=self.topofile_path)
        self.topofile_button.grid(column=0, row=0)
        
        # morsefile selection button
        self.morsefile = tk.Entry(self.inputs_frame, width=int(1.15*maxwidth), font=font_settings)
        self.morsefile.insert(0, morsefile)
        self.morsefile.grid(column=1, row=1)
        self.morsefile_button = tk.Button(self.inputs_frame, text='morsefile', font=font_settings, command=self.morsefile_path)
        self.morsefile_button.grid(column=0, row=1)

        # parent_directory entry
        self.parent_directory = tk.Entry(self.inputs_frame, width=int(1.15*maxwidth), font=font_settings)
        self.parent_directory.insert(0, parent_directory)
        self.parent_directory.grid(column=1, row=2)
        self.dir_button = tk.Button(self.inputs_frame, text='parent_directory', font=font_settings, command=self.directory_path)
        self.dir_button.grid(column=0, row=2)
        
        # newfile entry
        self.newfile = tk.Entry(self.inputs_frame, width=int(1.15*maxwidth), font=font_settings)
        self.newfile.insert(0, newfile)
        self.newfile.grid(column=1, row=3)
        self.newfile_label = tk.Label(self.inputs_frame, text='newfile', font=font_settings)
        self.newfile_label.grid(column=0, row=3)
        
        # ff_class drop down menu
        styles = ['1', '2']
        self.ff_class = ttk.Combobox(self.inputs_frame, values=styles, width=int(1.118*maxwidth), font=font_settings)
        self.ff_class.current(styles.index(ff_class))
        self.ff_class.grid(column=1, row=5)
        self.ff_class_label = tk.Label(self.inputs_frame, text='ff_class', font=font_settings)
        self.ff_class_label.grid(column=0, row=5)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.inputs_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=ypadding)
            
        
        #---------------#
        # Options frame #
        #---------------#
        # Initalize  options frame
        self.options_frame = tk.LabelFrame(self.frame2, text='Options', font=font_settings)
        self.options_frame.grid(row=1, column=0, columnspan=2, sticky='news', padx=xpadding, pady=ypadding)
        
        # atom_style drop down
        styles = ['full', 'charge', 'molecular', 'angle', 'bond', 'atomic', 'dipole', 'dpd', 'line']
        self.atom_style = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/7), font=font_settings)
        self.atom_style.current(styles.index(atom_style))
        self.atom_style.grid(column=0, row=1)
        self.atom_style_label = tk.Label(self.options_frame, text='atom_style', font=font_settings)
        self.atom_style_label.grid(column=0, row=0)
                
        # include_type_labels drop down menu
        styles = [True, False]
        self.include_type_labels = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/7), font=font_settings)
        self.include_type_labels.current(styles.index(include_type_labels))
        self.include_type_labels.grid(column=1, row=1)
        self.include_type_labels_label = tk.Label(self.options_frame, text='include_type_labels', font=font_settings)
        self.include_type_labels_label.grid(column=1, row=0)
                
        # include_rcut drop down menu
        styles = [True, False]
        self.include_rcut = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/7), font=font_settings)
        self.include_rcut.current(styles.index(include_rcut))
        self.include_rcut.grid(column=2, row=1)
        self.include_rcut_label = tk.Label(self.options_frame, text='include_rcut', font=font_settings)
        self.include_rcut_label.grid(column=2, row=0)
        
        # min_bond_length entry
        self.min_bond_length = tk.Entry(self.options_frame, width=int(maxwidth/5.75), font=font_settings)
        self.min_bond_length.insert(0, min_bond_length)
        self.min_bond_length.grid(column=3, row=1)
        self.min_bond_length_label = tk.Label(self.options_frame, text='min_bond_length', font=font_settings)
        self.min_bond_length_label.grid(column=3, row=0)
        
        # bondbreak_scale entry
        self.bondbreak_scale = tk.Entry(self.options_frame, width=int(maxwidth/5.75), font=font_settings)
        self.bondbreak_scale.insert(0, bondbreak_scale)
        self.bondbreak_scale.grid(column=4, row=1)
        self.bondbreak_scale_label = tk.Label(self.options_frame, text='bondbreak_scale', font=font_settings)
        self.bondbreak_scale_label.grid(column=4, row=0)
        
        # alpha_scale entry
        self.alpha_scale = tk.Entry(self.options_frame, width=int(maxwidth/5.75), font=font_settings)
        self.alpha_scale.insert(0, alpha_scale)
        self.alpha_scale.grid(column=5, row=1)
        self.alpha_scale_label = tk.Label(self.options_frame, text='alpha_scale', font=font_settings)
        self.alpha_scale_label.grid(column=5, row=0)
        
        # coeffs2skip entry
        self.coeffs2skip = tk.Entry(self.options_frame, width=int(0.9*maxwidth), font=font_settings)
        self.coeffs2skip.insert(0, ','.join([str(i) for i in coeffs2skip]))
        self.coeffs2skip.grid(column=2, row=2, columnspan=4)
        self.coeffs2skip_label = tk.Label(self.options_frame, text='Bond CoeffIDs 2 skip\n(comma separated w/no whitespace)', font=font_settings)
        self.coeffs2skip_label.grid(column=0, row=2, columnspan=2)

        # Add padding to all frames in self.inputs_frame
        for widget in self.options_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))
        
        #-------------------#
        # plt Options frame #
        #-------------------#
        # Initalize  plt_options frame
        self.plt_options_frame = tk.LabelFrame(self.frame2, text='Plot/Alphas2check Options', font=font_settings)
        self.plt_options_frame.grid(row=2, column=0, columnspan=2, sticky='news', padx=xpadding, pady=ypadding)
        
        # Radius spec start
        self.rss = tk.Entry(self.plt_options_frame, width=int(maxwidth/5.75), font=font_settings)
        self.rss.insert(0, radius_specs['start'])
        self.rss.grid(column=0, row=1)
        self.rss_label = tk.Label(self.plt_options_frame, text='plot start radius', font=font_settings)
        self.rss_label.grid(column=0, row=0)
        
        # Radius spec end
        self.rse = tk.Entry(self.plt_options_frame, width=int(maxwidth/5.75), font=font_settings)
        self.rse.insert(0, radius_specs['end'])
        self.rse.grid(column=1, row=1)
        self.rse_label = tk.Label(self.plt_options_frame, text='plot end radius', font=font_settings)
        self.rse_label.grid(column=1, row=0)
        
        # Radius spec inc
        self.rsi = tk.Entry(self.plt_options_frame, width=int(maxwidth/5.75), font=font_settings)
        self.rsi.insert(0, radius_specs['increment'])
        self.rsi.grid(column=2, row=1)
        self.rsi_label = tk.Label(self.plt_options_frame, text='radius increment', font=font_settings)
        self.rsi_label.grid(column=2, row=0)
        
        # Alpha spec start
        self.ass = tk.Entry(self.plt_options_frame, width=int(maxwidth/5.75), font=font_settings)
        self.ass.insert(0, alpha_specs['start'])
        self.ass.grid(column=3, row=1)
        self.ass_label = tk.Label(self.plt_options_frame, text='smallest alpha', font=font_settings)
        self.ass_label.grid(column=3, row=0)
        
        # Alpha spec end
        self.ase = tk.Entry(self.plt_options_frame, width=int(maxwidth/5.75), font=font_settings)
        self.ase.insert(0, alpha_specs['end'])
        self.ase.grid(column=4, row=1)
        self.ase_label = tk.Label(self.plt_options_frame, text='largest alpha', font=font_settings)
        self.ase_label.grid(column=4, row=0)
        
        # Alpha spec inc
        self.asi = tk.Entry(self.plt_options_frame, width=int(maxwidth/5.75), font=font_settings)
        self.asi.insert(0, alpha_specs['increment'])
        self.asi.grid(column=5, row=1)
        self.asi_label = tk.Label(self.plt_options_frame, text='alpha increment', font=font_settings)
        self.asi_label.grid(column=5, row=0)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.plt_options_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))
            
            
        #--------------------------------#
        # class2 crossterm Options frame #
        #--------------------------------#
        # Initalize  plt_options frame
        self.xterms_options_frame = tk.LabelFrame(self.frame1, text='Class2 crossterm Options', font=font_settings)
        self.xterms_options_frame.grid(row=3, column=0, columnspan=1, sticky='news', padx=xpadding, pady=ypadding)
        
        # zero_effected_xterms drop down menu
        styles = [True, False]
        self.zero_effected_xterms = ttk.Combobox(self.xterms_options_frame, values=styles, width=int(maxwidth/8.5), font=font_settings)
        self.zero_effected_xterms.current(styles.index(zero_effected_xterms))
        self.zero_effected_xterms.grid(column=0, row=1)
        self.zero_effected_xterms_label = tk.Label(self.xterms_options_frame, text='zero_effected_xterms', font=font_settings)
        self.zero_effected_xterms_label.grid(column=0, row=0)
        
        # class2xe_update drop down menu
        styles = [True, False]
        self.class2xe_update = ttk.Combobox(self.xterms_options_frame, values=styles, width=int(maxwidth/8.5), font=font_settings)
        self.class2xe_update.current(styles.index(class2xe_update))
        self.class2xe_update.grid(column=1, row=1)
        self.class2xe_update_label = tk.Label(self.xterms_options_frame, text='class2xe_update', font=font_settings)
        self.class2xe_update_label.grid(column=1, row=0)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.xterms_options_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))
            
            
        #--------------------#
        # file Options frame #
        #--------------------#
        # Initalize  plt_options frame
        self.file_options_frame = tk.LabelFrame(self.frame1, text='Files2write Options', font=font_settings)
        self.file_options_frame.grid(row=3, column=1, columnspan=1, sticky='news', padx=xpadding, pady=ypadding)
        
        # 'write_datafile' entry
        styles = [True, False]
        self.data = ttk.Combobox(self.file_options_frame, values=styles, width=int(maxwidth/8.5), font=font_settings)
        self.data.current(styles.index(files2write['write_datafile']))
        self.data.grid(column=0, row=1)
        self.data_label = tk.Label(self.file_options_frame, text='write_datafile', font=font_settings)
        self.data_label.grid(column=0, row=0)
        
        # 'write_pdffile' entry
        styles = [True, False]
        self.pdf = ttk.Combobox(self.file_options_frame, values=styles, width=int(maxwidth/8.5), font=font_settings)
        self.pdf.current(styles.index(files2write['write_pdffile']))
        self.pdf.grid(column=1, row=1)
        self.pdf_label = tk.Label(self.file_options_frame, text='write_pdffile', font=font_settings)
        self.pdf_label.grid(column=1, row=0)
        
        # 'write_bondbreak' entry
        styles = [True, False]
        self.bondbreak = ttk.Combobox(self.file_options_frame, values=styles, width=int(maxwidth/8.5), font=font_settings)
        self.bondbreak.current(styles.index(files2write['write_bondbreak']))
        self.bondbreak.grid(column=2, row=1)
        self.bondbreak_label = tk.Label(self.file_options_frame, text='write_bondbreak', font=font_settings)
        self.bondbreak_label.grid(column=2, row=0)
        
        # 'write_forcefield' entry
        styles = [True, False]
        self.forcefield = ttk.Combobox(self.file_options_frame, values=styles, width=int(maxwidth/8.5), font=font_settings)
        self.forcefield.current(styles.index(files2write['write_forcefield']))
        self.forcefield.grid(column=3, row=1)
        self.forcefield_label = tk.Label(self.file_options_frame, text='write_forcefield', font=font_settings)
        self.forcefield_label.grid(column=3, row=0)
        
        # Add padding to all frames in self.file_options_frame
        for widget in self.file_options_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))


        #------------#
        # Run button #
        #------------#
        self.run1 = tk.Button(self.frame1, text='Run LUNAR/auto_morse_bond_update.py', font=font_settings, command=self.run_LUNAR)
        self.run1.grid(row=4, column=0, columnspan=2, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        self.run2 = tk.Button(self.frame2, text='Run LUNAR/auto_morse_bond_update.py', font=font_settings, command=self.run_LUNAR)
        self.run2.grid(row=4, column=0, columnspan=2, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
        #-----------------#
        # update defaults #
        #-----------------#
        self.update1 = tk.Button(self.frame1, text='Save the current GUI settings as the default GUI settings', font=font_settings, command=self.update_py_script)
        self.update1.grid(row=5, column=0, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        self.update2 = tk.Button(self.frame2, text='Save the current GUI settings as the default GUI settings', font=font_settings, command=self.update_py_script)
        self.update2.grid(row=5, column=0, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
        #------------#
        # Quick help #
        #------------#
        self.quick_help1 = tk.Button(self.frame1, text='Quick help', font=font_settings, command=self.quickhelp)
        self.quick_help1.grid(row=5, column=1, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        self.quick_help2 = tk.Button(self.frame2, text='Quick help', font=font_settings, command=self.quickhelp)
        self.quick_help2.grid(row=5, column=1, sticky='news', padx=int(xpadding/2), pady=int(ypadding/2))
        
        
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
            txt = os.path.join(self.pwd, 'src/GUI_quick_help_pages/auto_morse_bond_update.txt')
            logged = []
            with open(txt, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    # Strip comment's and split by whitespace
                    line = line.split('#')[0]
                    line = line.rstrip()
                    logged.append(line)
        except: # except something failed
            logged.append('FAILED to read LUNAR/src/GUI_quick_help_pages/auto_morse_bond_update.txt document.')
            logged.append('Most likely cause is the auto_morse_bond_update.txt file was renamed in the')
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
    
    # Function to get filepath for morsefile
    def morsefile_path(self):
        ftypes = (('txt files', '*.txt'), ('all files', '*.*'))
        path = filedialog.askopenfilename(title='Open morsefile?', filetypes=ftypes)
        if path:
            path = os.path.relpath(path)
            self.morsefile.delete(0, tk.END); self.morsefile.insert(0, path);
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
        print('Terminating auto_morse_bond GUI'); self.root.destroy();
        return
    
    # Function to run imported LUNAR functions
    def run_LUNAR(self):  
        # Set up log
        log = io_functions.LUNAR_logger()
        log.configure(level='production', print2console=True,  write2log=True)
        valid_inputs = True
        
        # Get information from GUI
        boolean = {'False':False, 'True':True}
        topofile = self.topofile.get()
        morsefile = self.morsefile.get()
        parent_directory = self.parent_directory.get() 
        newfile = self.newfile.get()
        atom_style = self.atom_style.get()
        ff_class = self.ff_class.get()
        mass_map = self.mass_map
        include_type_labels = boolean[self.include_type_labels.get()]
        zero_effected_xterms = boolean[self.zero_effected_xterms.get()]
        try: radius_specs = {'start':float(self.rss.get()), 'end': float(self.rse.get()), 'increment':float(self.rsi.get())} 
        except: 
            log.GUI_error('ERROR a radius spec of start, end, or increment is not a float value.')
            valid_inputs = False
        try: alpha_specs  = {'start':float(self.ass.get()), 'end': float(self.ase.get()), 'increment':float(self.asi.get())} 
        except:
            log.GUI_error('ERROR an alpha spec of start, end, or increment is not a float value.')
            valid_inputs = False
        try: coeffs2skip = [int(i) for i in self.coeffs2skip.get().split(',')]
        except: coeffs2skip = []
        try: alpha_scale = float(self.alpha_scale.get())
        except:
            log.GUI_error('ERROR alpha_scale is not a float value.')
            valid_inputs = False
        try: bondbreak_scale = float(self.bondbreak_scale.get())
        except:
            log.GUI_error('ERROR bondbreak_scale is not a float value.')
            valid_inputs = False
        try: min_bond_length = float(self.min_bond_length.get())
        except:
            log.GUI_error('ERROR min_bond_length is not a float value.')
            valid_inputs = False
        files2write = {'write_datafile' : boolean[self.data.get()],
                       'write_pdffile'  : boolean[self.pdf.get()],
                       'write_bondbreak': boolean[self.bondbreak.get()],
                       'write_forcefield': boolean[self.forcefield.get()],
                       }
        include_rcut = boolean[self.include_rcut.get()]
        class2xe_update = boolean[self.class2xe_update.get()]

        # Run LUNAR/auto_morse_bond
        valid_inputs = True
        if valid_inputs:
            try: 
                inputs = (topofile, morsefile, parent_directory, newfile, mass_map, min_bond_length, coeffs2skip,
                          radius_specs, alpha_specs, alpha_scale, files2write, atom_style, zero_effected_xterms,
                          bondbreak_scale, ff_class, include_type_labels, class2xe_update, include_rcut, [], log)
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
        morsefile = io_functions.path_to_string(self.morsefile.get())
        newfile = self.newfile.get()
        atom_style = self.atom_style.get()
        ff_class = self.ff_class.get()
        include_type_labels = boolean[self.include_type_labels.get()]
        zero_effected_xterms = boolean[self.zero_effected_xterms.get()]
        radius_specs = {'start':float(self.rss.get()), 'end': float(self.rse.get()), 'increment':float(self.rsi.get())} 
        alpha_specs  = {'start':float(self.ass.get()), 'end': float(self.ase.get()), 'increment':float(self.asi.get())} 
        try: coeffs2skip = [int(i) for i in self.coeffs2skip.get().split(',')]
        except: coeffs2skip = []
        alpha_scale = float(self.alpha_scale.get())
        bondbreak_scale = float(self.bondbreak_scale.get())
        min_bond_length = float(self.min_bond_length.get())
        include_rcut = boolean[self.include_rcut.get()]
        files2write = {'write_datafile' : boolean[self.data.get()],
                       'write_pdffile'  : boolean[self.pdf.get()],
                       'write_bondbreak': boolean[self.bondbreak.get()],
                       'write_forcefield': boolean[self.forcefield.get()],
                       }                  
        class2xe_update = boolean[self.class2xe_update.get()]
        
        # Read current py script and re-write with new settings
        print('Updating settings in: {}, from current GUI settings'.format(self.filename))
        lines = psm.read(self.filename)
        with open(self.filename, 'w') as f:
            inputsflag = True; filesdict_flag = False
            for line in lines:
                # if line.startswith('use_GUI') and inputsflag:
                #     line = psm.parse_and_modify(line, True, stringflag=False, splitchar='=')
                # if line.startswith('topofile') and inputsflag:
                #     line = psm.parse_and_modify(line, topofile, stringflag=True, splitchar='=')
                if line.startswith('morsefile') and inputsflag:
                    line = psm.parse_and_modify(line, morsefile, stringflag=True, splitchar='=')
                if line.startswith('parent_directory') and inputsflag:
                    line = psm.parse_and_modify(line, parent_directory, stringflag=True, splitchar='=')
                if line.startswith('newfile') and inputsflag:
                    line = psm.parse_and_modify(line, newfile, stringflag=True, splitchar='=')
                if line.startswith('atom_style') and inputsflag:
                    line = psm.parse_and_modify(line, atom_style, stringflag=True, splitchar='=')
                if line.startswith('ff_class') and inputsflag:
                    line = psm.parse_and_modify(line, ff_class, stringflag=True, splitchar='=')
                if line.startswith('include_type_labels') and inputsflag:
                    line = psm.parse_and_modify(line, include_type_labels, stringflag=False, splitchar='=')
                if line.startswith('zero_effected_xterms') and inputsflag:
                    line = psm.parse_and_modify(line, zero_effected_xterms, stringflag=False, splitchar='=')
                if line.startswith('class2xe_update') and inputsflag:
                    line = psm.parse_and_modify(line, class2xe_update, stringflag=False, splitchar='=')
                if line.startswith('radius_specs') and inputsflag:
                    line = psm.parse_and_modify(line, str(radius_specs), stringflag=False, splitchar='=')
                if line.startswith('alpha_specs') and inputsflag:
                    line = psm.parse_and_modify(line, str(alpha_specs), stringflag=False, splitchar='=')
                if line.startswith('coeffs2skip') and inputsflag:
                    line = psm.parse_and_modify(line, str(coeffs2skip), stringflag=False, splitchar='=')
                if line.startswith('alpha_scale') and inputsflag:
                    line = psm.parse_and_modify(line, alpha_scale, stringflag=False, splitchar='=')
                if line.startswith('min_bond_length') and inputsflag:
                    line = psm.parse_and_modify(line, min_bond_length, stringflag=False, splitchar='=')
                if line.startswith('bondbreak_scale') and inputsflag:
                    line = psm.parse_and_modify(line, bondbreak_scale, stringflag=False, splitchar='=')
                if line.startswith('include_rcut') and inputsflag:
                    line = psm.parse_and_modify(line, include_rcut, stringflag=False, splitchar='=')   
                if line.startswith('files2write') and inputsflag: filesdict_flag = True
                if 'write_datafile' in line and inputsflag and not line.startswith('#') and filesdict_flag:
                    line = psm.parse_and_modify(line, files2write['write_datafile'], stringflag=False, splitchar=':')
                if 'write_pdffile' in line and inputsflag and not line.startswith('#') and filesdict_flag:
                    line = psm.parse_and_modify(line, files2write['write_pdffile'], stringflag=False, splitchar=':')
                if 'write_bondbreak' in line and inputsflag and not line.startswith('#') and filesdict_flag:
                    line = psm.parse_and_modify(line, files2write['write_bondbreak'], stringflag=False, splitchar=':')
                if 'write_forcefield' in line and inputsflag and not line.startswith('#') and filesdict_flag:
                    line = psm.parse_and_modify(line, files2write['write_forcefield'], stringflag=False, splitchar=':')        
                if line.startswith('if __name__ == "__main__":'): inputsflag = False
                f.write(line)
        return