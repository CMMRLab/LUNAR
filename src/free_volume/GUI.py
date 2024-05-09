# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 1st, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import src.io_functions as io_functions
from src.free_volume.main import main
import src.py_script_modifier as psm
from tkinter.scrolledtext import ScrolledText
from tkinter import filedialog
from tkinter import Toplevel
from tkinter import ttk
import tkinter as tk
import math
import os


#########################
# LUNAR/free_volume GUI #
#########################
class free_volume_GUI:
    def __init__(self, topofile, max_voxel_size, mass_map, vdw_radius, boundary, parent_directory,
                 compute_free_volume_distributions, files2write, run_mode, probe_diameter, vdw_method,
                 CUDA_threads_per_block_atoms, CUDA_threads_per_block_voxels, GUI_zoom):
        
        # Pass certain inputs as attribute (These will only be able to be adjusted from the python script)
        self.mass_map = mass_map
        self.vdw_radius = vdw_radius
        
        # Find present working directory
        self.pwd = os.getcwd()
        self.filepath = self.pwd
        self.filename = os.path.join(self.pwd, 'free_volume.py')

        # Initialize window
        self.root = tk.Tk()
        self.root.title('LUNAR/free_volume.py GUI v1.0')
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
        self.maxwidth = 90
        
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

        # parent_directory entry
        self.parent_directory = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.parent_directory.insert(0, parent_directory)
        self.parent_directory.grid(column=1, row=1)
        self.dir_button = tk.Button(self.inputs_frame, text='parent_directory', font=font_settings, command=self.directory_path)
        self.dir_button.grid(column=0, row=1)
        
        # run mode drop menu
        styles = ['stl', 'numpy', 'numba', 'numba-p', 'stl-dd', 'numba-dd', 'numba-ddp', 'CUDA', 'CUDA-dd']
        self.run_mode = ttk.Combobox(self.inputs_frame, values=styles, width=int(maxwidth/1.023), font=font_settings)
        self.run_mode.current(styles.index(run_mode))
        self.run_mode.grid(column=1, row=2)
        self.run_mode_label = tk.Label(self.inputs_frame, text='run_mode', font=font_settings)
        self.run_mode_label.grid(column=0, row=2)
        
        # CUDA_threads_per_block_atoms down menu
        m8 = [8, 16, 32, 64, 128, 256, 512, 1024]
        self.CUDA_threads_per_block_atoms = ttk.Combobox(self.inputs_frame, values=m8, width=int(maxwidth/1.023), font=font_settings)
        self.CUDA_threads_per_block_atoms.current(m8.index(CUDA_threads_per_block_atoms))
        self.CUDA_threads_per_block_atoms.grid(column=1, row=3)
        self.CUDA_threads_per_block_atoms_label = tk.Label(self.inputs_frame, text='CUDA_threads_per_block_atoms', font=font_settings)
        self.CUDA_threads_per_block_atoms_label.grid(column=0, row=3)
        
        # CUDA_threads_per_block_voxels down menu
        m8.append(0)
        self.CUDA_threads_per_block_voxels = ttk.Combobox(self.inputs_frame, values=m8, width=int(maxwidth/1.023), font=font_settings)
        self.CUDA_threads_per_block_voxels.current(m8.index(CUDA_threads_per_block_voxels))
        self.CUDA_threads_per_block_voxels.grid(column=1, row=4)
        self.CUDA_threads_per_block_voxels_label = tk.Label(self.inputs_frame, text='CUDA_threads_per_block_voxels', font=font_settings)
        self.CUDA_threads_per_block_voxels_label.grid(column=0, row=4)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.inputs_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=ypadding)
            
        
        #---------------#
        # Options frame #
        #---------------#
        # Initalize  options frame
        self.options_frame = tk.LabelFrame(self.frame, text='Options', font=font_settings)
        self.options_frame.grid(row=1, column=0, sticky='news', padx=xpadding, pady=ypadding)
        
        # boundary entry
        self.boundary = tk.Entry(self.options_frame, width=int(maxwidth/9), font=font_settings)
        self.boundary.insert(0, boundary)
        self.boundary.grid(column=0, row=1)
        self.boundary_label = tk.Label(self.options_frame, text='boundary', font=font_settings)
        self.boundary_label.grid(column=0, row=0)
        
        # max_voxel_size entry
        self.max_voxel_size = tk.Entry(self.options_frame, width=int(maxwidth/9), font=font_settings)
        self.max_voxel_size.insert(0, max_voxel_size)
        self.max_voxel_size.grid(column=1, row=1)
        self.max_voxel_size_label = tk.Label(self.options_frame, text='max_voxel_size', font=font_settings)
        self.max_voxel_size_label.grid(column=1, row=0)
        
        # probe_diameter entry
        self.probe_diameter = tk.Entry(self.options_frame, width=int(maxwidth/9), font=font_settings)
        self.probe_diameter.insert(0, probe_diameter)
        self.probe_diameter.grid(column=2, row=1)
        self.probe_diameter_label = tk.Label(self.options_frame, text='probe_diameter', font=font_settings)
        self.probe_diameter_label.grid(column=2, row=0)
        
        # compute_free_volume_distributions drop down menu
        styles = ['dict', 'class1', 'class2']
        self.vdw_method = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/4), font=font_settings)
        self.vdw_method.current(styles.index(vdw_method))
        self.vdw_method.grid(column=3, row=1)
        self.vdw_method_label = tk.Label(self.options_frame, text='vdw_method', font=font_settings)
        self.vdw_method_label.grid(column=3, row=0)
        
        # compute_free_volume_distributions drop down menu
        styles = [True, False]
        self.compute_free_volume_distributions = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/4), font=font_settings)
        self.compute_free_volume_distributions.current(styles.index(compute_free_volume_distributions))
        self.compute_free_volume_distributions.grid(column=4, row=1)
        self.compute_free_volume_distributions_label = tk.Label(self.options_frame, text='compute_free_volume_distributions', font=font_settings)
        self.compute_free_volume_distributions_label.grid(column=4, row=0)

        
        # Add padding to all frames in self.inputs_frame
        for widget in self.options_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))
            
        #--------------------#
        # file Options frame #
        #--------------------#
        # Initalize  plt_options frame
        self.file_options_frame = tk.LabelFrame(self.frame, text='Files2write Options', font=font_settings)
        self.file_options_frame.grid(row=3, column=0, sticky='news', padx=xpadding, pady=ypadding)
        
        # 'write_atoms_free' entry
        styles = [True, False]
        self.atoms_free = ttk.Combobox(self.file_options_frame, values=styles, width=int(maxwidth/4), font=font_settings)
        self.atoms_free.current(styles.index(files2write['write_atoms_free']))
        self.atoms_free.grid(column=0, row=1)
        self.atoms_free_label = tk.Label(self.file_options_frame, text='write_atoms_free', font=font_settings)
        self.atoms_free_label.grid(column=0, row=0)
        
        # 'write_bonds_free' entry
        styles = [True, False]
        self.bonds_free = ttk.Combobox(self.file_options_frame, values=styles, width=int(maxwidth/4), font=font_settings)
        self.bonds_free.current(styles.index(files2write['write_bonds_free']))
        self.bonds_free.grid(column=1, row=1)
        self.bonds_free_label = tk.Label(self.file_options_frame, text='write_bonds_free', font=font_settings)
        self.bonds_free_label.grid(column=1, row=0)
        
        # 'write_atoms_only' entry
        styles = [True, False]
        self.atoms_only = ttk.Combobox(self.file_options_frame, values=styles, width=int(maxwidth/4), font=font_settings)
        self.atoms_only.current(styles.index(files2write['write_atoms_only']))
        self.atoms_only.grid(column=2, row=1)
        self.atoms_only_label = tk.Label(self.file_options_frame, text='write_atoms_only', font=font_settings)
        self.atoms_only_label.grid(column=2, row=0)
        
        # 'write_free_only' entry
        styles = [True, False]
        self.free_only = ttk.Combobox(self.file_options_frame, values=styles, width=int(maxwidth/4), font=font_settings)
        self.free_only.current(styles.index(files2write['write_free_only' ]))
        self.free_only.grid(column=3, row=1)
        self.free_only_label = tk.Label(self.file_options_frame, text='write_free_only', font=font_settings)
        self.free_only_label.grid(column=3, row=0)
        
        # 'write_spat_dis-x' entry
        styles = [True, False]
        self.spat_dis_x = ttk.Combobox(self.file_options_frame, values=styles, width=int(maxwidth/4), font=font_settings)
        self.spat_dis_x.current(styles.index(files2write['write_spat_dis-x']))
        self.spat_dis_x.grid(column=0, row=3)
        self.spat_dis_x_label = tk.Label(self.file_options_frame, text='write_spat_dis-x', font=font_settings)
        self.spat_dis_x_label.grid(column=0, row=2)
        
        # 'write_spat_dis-y' entry
        styles = [True, False]
        self.spat_dis_y = ttk.Combobox(self.file_options_frame, values=styles, width=int(maxwidth/4), font=font_settings)
        self.spat_dis_y.current(styles.index(files2write['write_spat_dis-y']))
        self.spat_dis_y.grid(column=1, row=3)
        self.spat_dis_y_label = tk.Label(self.file_options_frame, text='write_spat_dis-y', font=font_settings)
        self.spat_dis_y_label.grid(column=1, row=2)
        
        # 'write_spat_dis-z' entry
        styles = [True, False]
        self.spat_dis_z = ttk.Combobox(self.file_options_frame, values=styles, width=int(maxwidth/4), font=font_settings)
        self.spat_dis_z.current(styles.index(files2write['write_spat_dis-z']))
        self.spat_dis_z.grid(column=2, row=3)
        self.spat_dis_z_label = tk.Label(self.file_options_frame, text='write_spat_dis-z', font=font_settings)
        self.spat_dis_z_label.grid(column=2, row=2)
        
        # 'write_all_voxels' entry
        styles = [True, False]
        self.all_voxels = ttk.Combobox(self.file_options_frame, values=styles, width=int(maxwidth/4), font=font_settings)
        self.all_voxels.current(styles.index(files2write['write_all_voxels']))
        self.all_voxels.grid(column=3, row=3)
        self.all_voxels_label = tk.Label(self.file_options_frame, text='write_all_voxels', font=font_settings)
        self.all_voxels_label.grid(column=3, row=2)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.file_options_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))

        #------------#
        # Run button #
        #------------#
        self.run = tk.Button(self.frame, text='Run LUNAR/free_volume.py', font=font_settings, command=self.run_LUNAR)
        self.run.grid(row=4, column=0, sticky='news', padx=xpadding, pady=ypadding)
        
        #-----------------#
        # update defaults #
        #-----------------#
        self.update = tk.Button(self.frame, text='Save the current GUI settings as the default GUI settings', font=font_settings, command=self.update_py_script)
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
    
    # Function to get directory
    def directory_path(self):
        path =filedialog.askdirectory(initialdir=self.pwd)
        if path:
            path = os.path.relpath(path)
            self.parent_directory.delete(0, tk.END); self.parent_directory.insert(0, path);
        return
    
    # Closing command    
    def closing(self):
        print('Terminating free_volume GUI'); self.root.destroy();
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
        parent_directory = self.parent_directory.get() 
        mass_map = self.mass_map
        vdw_radius = self.vdw_radius
        boundary = self.boundary.get()
        try: max_voxel_size = float(self.max_voxel_size.get())
        except:
            valid_inputs = False
            log.GUI_error('ERROR max_voxel_size is not float valued')
        try: probe_diameter = float(self.probe_diameter.get())
        except: 
            valid_inputs = False
            log.GUI_error('ERROR probe_diameter is not float valued')
        compute_free_volume_distributions = boolean[self.compute_free_volume_distributions.get()]
        files2write = {'write_atoms_free' : boolean[self.atoms_free.get()],
                       'write_bonds_free' : boolean[self.bonds_free.get()],
                       'write_atoms_only' : boolean[self.atoms_only.get()],
                       'write_free_only'  : boolean[self.free_only.get()],
                       'write_all_voxels' : boolean[self.all_voxels.get()],
                       'write_spat_dis-x' : boolean[self.spat_dis_x.get()],
                       'write_spat_dis-y' : boolean[self.spat_dis_y.get()],
                       'write_spat_dis-z' : boolean[self.spat_dis_z.get()],
                       }
        run_mode = self.run_mode.get()
        vdw_method = self.vdw_method.get()
        CUDA_threads_per_block_atoms = int(self.CUDA_threads_per_block_atoms.get())
        CUDA_threads_per_block_voxels = int(self.CUDA_threads_per_block_voxels.get())
        
        # Run LUNAR/free_volume
        if valid_inputs:
            try: main(topofile, max_voxel_size, mass_map, vdw_radius, boundary, parent_directory,
                      compute_free_volume_distributions, files2write, run_mode, probe_diameter,
                      vdw_method, CUDA_threads_per_block_atoms, CUDA_threads_per_block_voxels, [], log=log)
            except Exception as error:
                log.GUI_error('{}: {}'.format(type(error).__name__, str(error)))
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
        # Get information from GUI
        boolean = {'False':False, 'True':True}
        topofile = io_functions.path_to_string(self.topofile.get())
        parent_directory = io_functions.path_to_string(self.parent_directory.get()) 
        boundary = self.boundary.get()
        max_voxel_size = float(self.max_voxel_size.get())
        try:
            probe_diameter = float(self.probe_diameter.get())
            probe_type = 'float'
        except: 
            probe_diameter = self.probe_diameter.get()
            probe_type = 'string'
        compute_free_volume_distributions = boolean[self.compute_free_volume_distributions.get()]
        files2write = {'write_atoms_free' : boolean[self.atoms_free.get()],
                       'write_bonds_free' : boolean[self.bonds_free.get()],
                       'write_atoms_only' : boolean[self.atoms_only.get()],
                       'write_free_only'  : boolean[self.free_only.get()],
                       'write_all_voxels' : boolean[self.all_voxels.get()],
                       'write_spat_dis-x' : boolean[self.spat_dis_x.get()],
                       'write_spat_dis-y' : boolean[self.spat_dis_y.get()],
                       'write_spat_dis-z' : boolean[self.spat_dis_z.get()],
                       }
        run_mode = self.run_mode.get()
        vdw_method = self.vdw_method.get()
        CUDA_threads_per_block_atoms = int(self.CUDA_threads_per_block_atoms.get())
        CUDA_threads_per_block_voxels = int(self.CUDA_threads_per_block_voxels.get())
        
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
                if line.startswith('parent_directory') and inputsflag:
                    line = psm.parse_and_modify(line, parent_directory, stringflag=True, splitchar='=')
                if line.startswith('boundary') and inputsflag:
                     line = psm.parse_and_modify(line, boundary, stringflag=True, splitchar='=')  
                if line.startswith('run_mode') and inputsflag:
                     line = psm.parse_and_modify(line, run_mode, stringflag=True, splitchar='=') 
                if line.startswith('CUDA_threads_per_block_atoms') and inputsflag:
                    line = psm.parse_and_modify(line, CUDA_threads_per_block_atoms, stringflag=False, splitchar='=')
                if line.startswith('CUDA_threads_per_block_voxels') and inputsflag:
                    line = psm.parse_and_modify(line, CUDA_threads_per_block_voxels, stringflag=False, splitchar='=')
                if line.startswith('max_voxel_size') and inputsflag:
                    line = psm.parse_and_modify(line, max_voxel_size, stringflag=False, splitchar='=')
                if line.startswith('vdw_method') and inputsflag:
                     line = psm.parse_and_modify(line, vdw_method, stringflag=True, splitchar='=')  
                if line.startswith('probe_diameter') and inputsflag:
                    if probe_type == 'float': stringflag=False
                    else: stringflag=True
                    line = psm.parse_and_modify(line, probe_diameter, stringflag=stringflag, splitchar='=')
                if line.startswith('compute_free_volume_distributions') and inputsflag:
                    line = psm.parse_and_modify(line, compute_free_volume_distributions, stringflag=False, splitchar='=')
                if line.startswith('files2write') and inputsflag: filesdict_flag = True
                
                if 'write_atoms_free' in line and inputsflag and not line.startswith('#') and filesdict_flag:
                    line = psm.parse_and_modify(line, files2write['write_atoms_free'], stringflag=False, splitchar=':')
                    
                if 'write_bonds_free' in line and inputsflag and not line.startswith('#') and filesdict_flag:
                    line = psm.parse_and_modify(line, files2write['write_bonds_free'], stringflag=False, splitchar=':')
                
                
                if 'write_atoms_only' in line and inputsflag and not line.startswith('#') and filesdict_flag:
                    line = psm.parse_and_modify(line, files2write['write_atoms_only'], stringflag=False, splitchar=':')
                if 'write_free_only' in line and inputsflag and not line.startswith('#') and filesdict_flag:
                    line = psm.parse_and_modify(line, files2write['write_free_only'], stringflag=False, splitchar=':')
                if 'write_all_voxels' in line and inputsflag and not line.startswith('#') and filesdict_flag:
                    line = psm.parse_and_modify(line, files2write['write_all_voxels'], stringflag=False, splitchar=':')
                if 'write_spat_dis-x' in line and inputsflag and not line.startswith('#') and filesdict_flag:
                    line = psm.parse_and_modify(line, files2write['write_spat_dis-x'], stringflag=False, splitchar=':')
                if 'write_spat_dis-y' in line and inputsflag and not line.startswith('#') and filesdict_flag:
                    line = psm.parse_and_modify(line, files2write['write_spat_dis-y'], stringflag=False, splitchar=':')
                if 'write_spat_dis-z' in line and inputsflag and not line.startswith('#') and filesdict_flag:
                    line = psm.parse_and_modify(line, files2write['write_spat_dis-z'], stringflag=False, splitchar=':')
                    
                    
                if line.startswith('if __name__ == "__main__":'): inputsflag = False
                f.write(line)
        return