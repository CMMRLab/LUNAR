# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 5th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
from src.auto_morse_bond.main import main
import src.io_functions as io_functions
import src.py_script_modifier as psm
from tkinter.scrolledtext import ScrolledText
from tkinter import filedialog
from tkinter import Toplevel
from tkinter import ttk
import tkinter as tk
import math
import os


#############################
# LUNAR/auto_morse_bond GUI #
#############################
class auto_morse_bond_GUI:
    def __init__(self, topofile, morsefile, parent_directory, newfile, mass_map, min_bond_length, coeffs2skip,
                 radius_specs, alpha_specs, alpha_scale, files2write, atom_style, zero_effected_xterms,
                 bondbreak_scale, ff_class, include_type_labels, include_rcut, GUI_zoom):
        
        # Pass certain inputs as attribute (These will only be able to be adjusted from the python script)
        self.mass_map = mass_map
        
        # Find present working directory
        self.pwd = os.getcwd()
        self.filepath = self.pwd
        self.filename = os.path.join(self.pwd, 'auto_morse_bond_update.py')

        # Initialize window
        self.root = tk.Tk()
        self.root.title('LUNAR/auto_morse_bond_update.py GUI v1.0')
        self.root.resizable(width=False, height=False)
        
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
        self.inputs_frame.grid(row=0, column=0, padx=xpadding, pady=ypadding)
        
        # topofile selection button
        self.topofile = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.topofile.insert(0, topofile)
        self.topofile.grid(column=1, row=0)
        self.topofile_button = tk.Button(self.inputs_frame, text='topofile\n(LAMMPS datafile to insert morse bonds)', font=font_settings, command=self.topofile_path)
        self.topofile_button.grid(column=0, row=0)
        
        # morsefile selection button
        self.morsefile = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.morsefile.insert(0, morsefile)
        self.morsefile.grid(column=1, row=1)
        self.morsefile_button = tk.Button(self.inputs_frame, text='morsefile\n(file containing bond typing rules\nand dissociation energy parameters)', font=font_settings, command=self.morsefile_path)
        self.morsefile_button.grid(column=0, row=1)

        # parent_directory entry
        self.parent_directory = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.parent_directory.insert(0, parent_directory)
        self.parent_directory.grid(column=1, row=2)
        self.dir_button = tk.Button(self.inputs_frame, text='parent_directory\n(location where files will be written to)', font=font_settings, command=self.directory_path)
        self.dir_button.grid(column=0, row=2)
        
        # newfile entry
        self.newfile = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.newfile.insert(0, newfile)
        self.newfile.grid(column=1, row=3)
        self.newfile_label = tk.Label(self.inputs_frame, text='newfile\n(appends to basename of topofile\nto set written files basename)', font=font_settings)
        self.newfile_label.grid(column=0, row=3)

        # ff_class drop down menu
        styles = [1, 2]
        self.ff_class = ttk.Combobox(self.inputs_frame, values=styles, width=maxwidth-3, font=font_settings)
        self.ff_class.current(styles.index(ff_class))
        self.ff_class.grid(column=1, row=4)
        self.ff_class_label = tk.Label(self.inputs_frame, text='ff_class\n(e.g. PCFF, COMPASS, CFF91 = 2\nCHARMM, OPLS, DREIDING, CVFF = 1)', font=font_settings)
        self.ff_class_label.grid(column=0, row=4)
        
        # bondbreak_scale entry
        self.bondbreak_scale = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.bondbreak_scale.insert(0, bondbreak_scale)
        self.bondbreak_scale.grid(column=1, row=5)
        self.bondbreak_scale_label = tk.Label(self.inputs_frame, text='bondbreak_scale\n(multiple of equilibrium bond length, default = 1.75)', font=font_settings)
        self.bondbreak_scale_label.grid(column=0, row=5)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.inputs_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=ypadding)
            
        
        #---------------#
        # Options frame #
        #---------------#
        # Initalize  options frame
        self.options_frame = tk.LabelFrame(self.frame2, text='Options', font=font_settings)
        self.options_frame.grid(row=1, column=0, sticky='news', padx=xpadding, pady=ypadding)
        
        # atom_style drop down
        styles = ['full', 'charge', 'molecular', 'angle', 'bond', 'atomic', 'dipole', 'dpd', 'line']
        self.atom_style = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/10), font=font_settings)
        self.atom_style.current(styles.index(atom_style))
        self.atom_style.grid(column=1, row=1)
        self.atom_style_label = tk.Label(self.options_frame, text='atom_style\n(formats Atoms section style)', font=font_settings)
        self.atom_style_label.grid(column=1, row=0)
        
        # zero_effected_xterms drop down menu
        styles = [True, False]
        self.zero_effected_xterms = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/10), font=font_settings)
        self.zero_effected_xterms.current(styles.index(zero_effected_xterms))
        self.zero_effected_xterms.grid(column=2, row=1)
        self.zero_effected_xterms_label = tk.Label(self.options_frame, text='zero_effected_xterms\n', font=font_settings)
        self.zero_effected_xterms_label.grid(column=2, row=0)
                
        # include_type_labels drop down menu
        styles = [True, False]
        self.include_type_labels = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/10), font=font_settings)
        self.include_type_labels.current(styles.index(include_type_labels))
        self.include_type_labels.grid(column=3, row=1)
        self.include_type_labels_label = tk.Label(self.options_frame, text='include_type_labels\n', font=font_settings)
        self.include_type_labels_label.grid(column=3, row=0)
                
        # include_rcut drop down menu
        styles = [True, False]
        self.include_rcut = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/10), font=font_settings)
        self.include_rcut.current(styles.index(include_rcut))
        self.include_rcut.grid(column=4, row=1)
        self.include_rcut_label = tk.Label(self.options_frame, text='include_rcut\n(adds rcut value to shift morse potential)', font=font_settings)
        self.include_rcut_label.grid(column=4, row=0)
        
        # min_bond_length entry
        self.min_bond_length = tk.Entry(self.options_frame, width=int(maxwidth/10), font=font_settings)
        self.min_bond_length.insert(0, min_bond_length)
        self.min_bond_length.grid(column=5, row=1)
        self.min_bond_length_label = tk.Label(self.options_frame, text='minimum bond r0\n to update to morse', font=font_settings)
        self.min_bond_length_label.grid(column=5, row=0)
        
        # alpha_scale entry
        self.alpha_scale = tk.Entry(self.options_frame, width=int(maxwidth/10), font=font_settings)
        self.alpha_scale.insert(0, alpha_scale)
        self.alpha_scale.grid(column=1, row=3)
        self.alpha_scale_label = tk.Label(self.options_frame, text='alpha_scale\n(multiple of best fit alpha, default = 1.0', font=font_settings)
        self.alpha_scale_label.grid(column=1, row=2)
        
        # coeffs2skip entry
        self.coeffs2skip = tk.Entry(self.options_frame, width=maxwidth, font=font_settings)
        self.coeffs2skip.insert(0, ','.join([str(i) for i in coeffs2skip]))
        self.coeffs2skip.grid(column=2, row=3, columnspan=4)
        self.coeffs2skip_label = tk.Label(self.options_frame, text='Bond CoeffIDs 2 skip\n(comma separated w/no whitespace)', font=font_settings)
        self.coeffs2skip_label.grid(column=2, row=2, columnspan=4)

        
        # Add padding to all frames in self.inputs_frame
        for widget in self.options_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))
        
        #-------------------#
        # plt Options frame #
        #-------------------#
        # Initalize  plt_options frame
        self.plt_options_frame = tk.LabelFrame(self.frame2, text='Plot/Alphas2check Options', font=font_settings)
        self.plt_options_frame.grid(row=2, column=0, sticky='news', padx=xpadding, pady=ypadding)
        
        # Radius spec start
        self.rss = tk.Entry(self.plt_options_frame, width=int(maxwidth/6.67), font=font_settings)
        self.rss.insert(0, radius_specs['start'])
        self.rss.grid(column=0, row=1)
        self.rss_label = tk.Label(self.plt_options_frame, text='plot start radius\n(default = 0.0)', font=font_settings)
        self.rss_label.grid(column=0, row=0)
        
        # Radius spec end
        self.rse = tk.Entry(self.plt_options_frame, width=int(maxwidth/6.67), font=font_settings)
        self.rse.insert(0, radius_specs['end'])
        self.rse.grid(column=1, row=1)
        self.rse_label = tk.Label(self.plt_options_frame, text='plot end radius\n(default = 8.0)', font=font_settings)
        self.rse_label.grid(column=1, row=0)
        
        # Radius spec inc
        self.rsi = tk.Entry(self.plt_options_frame, width=int(maxwidth/6.67), font=font_settings)
        self.rsi.insert(0, radius_specs['increment'])
        self.rsi.grid(column=2, row=1)
        self.rsi_label = tk.Label(self.plt_options_frame, text='radius increment\n(default = 0.01)', font=font_settings)
        self.rsi_label.grid(column=2, row=0)
        
        # Alpha spec start
        self.ass = tk.Entry(self.plt_options_frame, width=int(maxwidth/6.67), font=font_settings)
        self.ass.insert(0, alpha_specs['start'])
        self.ass.grid(column=3, row=1)
        self.ass_label = tk.Label(self.plt_options_frame, text='smallest possible alpha\n(default = 0.5)', font=font_settings)
        self.ass_label.grid(column=3, row=0)
        
        # Alpha spec end
        self.ase = tk.Entry(self.plt_options_frame, width=int(maxwidth/6.67), font=font_settings)
        self.ase.insert(0, alpha_specs['end'])
        self.ase.grid(column=4, row=1)
        self.ase_label = tk.Label(self.plt_options_frame, text='largest possible alpha\n(default = 3.5)', font=font_settings)
        self.ase_label.grid(column=4, row=0)
        
        # Alpha spec inc
        self.asi = tk.Entry(self.plt_options_frame, width=int(maxwidth/6.67), font=font_settings)
        self.asi.insert(0, alpha_specs['increment'])
        self.asi.grid(column=5, row=1)
        self.asi_label = tk.Label(self.plt_options_frame, text='alpha increment\n(default = 0.1)', font=font_settings)
        self.asi_label.grid(column=5, row=0)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.plt_options_frame.winfo_children():
            widget.grid_configure(padx=xpadding+2, pady=int(ypadding/2))
            
            
        #--------------------#
        # file Options frame #
        #--------------------#
        # Initalize  plt_options frame
        self.file_options_frame = tk.LabelFrame(self.frame2, text='Files2write Options', font=font_settings)
        self.file_options_frame.grid(row=3, column=0, sticky='news', padx=xpadding, pady=ypadding)
        
        # 'write_datafile' entry
        styles = [True, False]
        self.data = ttk.Combobox(self.file_options_frame, values=styles, width=int(maxwidth/9), font=font_settings)
        self.data.current(styles.index(files2write['write_datafile']))
        self.data.grid(column=0, row=1)
        self.data_label = tk.Label(self.file_options_frame, text='write_datafile', font=font_settings)
        self.data_label.grid(column=0, row=0)
        
        # 'write_pdffile' entry
        styles = [True, False]
        self.pdf = ttk.Combobox(self.file_options_frame, values=styles, width=int(maxwidth/9), font=font_settings)
        self.pdf.current(styles.index(files2write['write_pdffile']))
        self.pdf.grid(column=1, row=1)
        self.pdf_label = tk.Label(self.file_options_frame, text='write_pdffile', font=font_settings)
        self.pdf_label.grid(column=1, row=0)
        
        # 'write_bondbreak' entry
        styles = [True, False]
        self.bondbreak = ttk.Combobox(self.file_options_frame, values=styles, width=int(maxwidth/9), font=font_settings)
        self.bondbreak.current(styles.index(files2write['write_bondbreak']))
        self.bondbreak.grid(column=2, row=1)
        self.bondbreak_label = tk.Label(self.file_options_frame, text='write_bondbreak', font=font_settings)
        self.bondbreak_label.grid(column=2, row=0)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.file_options_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))


        #------------#
        # Run button #
        #------------#
        self.run = tk.Button(self.frame1, text='Run LUNAR/auto_morse_bond_update.py', font=font_settings, command=self.run_LUNAR)
        self.run.grid(row=4, column=0, sticky='news', padx=xpadding, pady=ypadding)
        
        #-----------------#
        # update defaults #
        #-----------------#
        self.update = tk.Button(self.frame1, text='Save the current GUI settings as the default GUI settings', font=font_settings, command=self.update_py_script)
        self.update.grid(row=5, column=0, sticky='news', padx=xpadding, pady=ypadding)
        
        
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
        ftypes = (('data files', '*.data *.data.gz'), ('all files', '*.*'))
        path = filedialog.askopenfilename(initialdir=self.filepath, title='Open topofile?', filetypes=ftypes)
        if path:
            self.filepath = os.path.dirname(os.path.abspath(path))
            path = os.path.relpath(path)
            self.topofile.delete(0, tk.END); self.topofile.insert(0, path);
        return
    
    # Function to get filepath for morsefile
    def morsefile_path(self):
        ftypes = (('txt files', '*.txt'), ('all files', '*.*'))
        path = filedialog.askopenfilename(initialdir=self.filepath, title='Open morsefile?', filetypes=ftypes)
        if path:
            self.filepath = os.path.dirname(os.path.abspath(path))
            path = os.path.relpath(path)
            self.morsefile.delete(0, tk.END); self.morsefile.insert(0, path);
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
        print('Terminating auto_morse_bond GUI'); self.root.destroy();
        return
    
    # Function to run imported LUNAR functions
    def run_LUNAR(self):  
        # Set up log
        log = io_functions.LUNAR_logger()
        log.configure(level='production', print2console=True,  write2log=True)
        valid_inputs = True
        
        # Get information from GUI
        tk2ff = {'0':0, '1':1, '2':2, 'r':'r', 'd':'d', 's1':'s1', 's2':'s2'}
        boolean = {'False':False, 'True':True}
        topofile = self.topofile.get()
        morsefile = self.morsefile.get()
        parent_directory = self.parent_directory.get() 
        newfile = self.newfile.get()
        atom_style = self.atom_style.get()
        ff_class = tk2ff[self.ff_class.get()]
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
                       'write_bondbreak': boolean[self.bondbreak.get()]}
        include_rcut = boolean[self.include_rcut.get()]


        # Run LUNAR/auto_morse_bond
        valid_inputs = True
        if valid_inputs:
            try: main(topofile, morsefile, parent_directory, newfile, mass_map, min_bond_length, coeffs2skip,
                      radius_specs, alpha_specs, alpha_scale, files2write, atom_style, zero_effected_xterms,
                      bondbreak_scale, ff_class, include_type_labels, include_rcut, log=log)
            except: pass
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
        tk2ff = {'0':0, '1':1, '2':2, 'r':'r', 'd':'d', 's1':'s1', 's2':'s2'}
        boolean = {'False':False, 'True':True}
        topofile = self.topofile.get()
        parent_directory = self.parent_directory.get() 
        morsefile = self.morsefile.get()
        newfile = self.newfile.get()
        atom_style = self.atom_style.get()
        ff_class = tk2ff[self.ff_class.get()]
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
                       'write_bondbreak': boolean[self.bondbreak.get()]}
        
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
                    line = psm.parse_and_modify(line, ff_class, stringflag=False, splitchar='=')
                if line.startswith('include_type_labels') and inputsflag:
                    line = psm.parse_and_modify(line, include_type_labels, stringflag=False, splitchar='=')
                if line.startswith('zero_effected_xterms') and inputsflag:
                    line = psm.parse_and_modify(line, zero_effected_xterms, stringflag=False, splitchar='=')
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
                if line.startswith('if __name__ == "__main__":'): inputsflag = False
                f.write(line)
        return