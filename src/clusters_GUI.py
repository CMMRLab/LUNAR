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
import src.io_functions as io_functions
import src.py_script_modifier as psm
import src.clusters as clusters
from tkinter.scrolledtext import ScrolledText
from tkinter import filedialog
from tkinter import Toplevel
from tkinter import ttk
import tkinter as tk
import threading
import traceback
import math
import glob
import time
import os




###################################
# LUNAR/auto_cluster_analysis GUI #
###################################
class GUI:
    def __init__(self, topofile, N0, txtfile, fav, GUI_zoom, pflag=True):
        self.pflag = pflag
        
        # Find present working directory
        self.pwd = os.getcwd()
        self.filepath = self.pwd
        self.filename = os.path.join(self.pwd, 'cluster_analysis.py')

        # Initialize window
        self.root = tk.Tk()
        self.root.title('LUNAR/cluster_analysis.py GUI v1.0')
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
        self.topofile = tk.Entry(self.inputs_frame, font=font_settings, width=maxwidth)
        self.topofile.insert(0, topofile)
        self.topofile.grid(column=1, row=0)
        self.topofile_button = tk.Button(self.inputs_frame, text='topofile', font=font_settings, command=self.topofile_path)
        self.topofile_button.grid(column=0, row=0)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.inputs_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=ypadding)
            
        
        #---------------#
        # Options frame #
        #---------------#
        # Initalize  options frame
        self.options_frame = tk.LabelFrame(self.frame, text='Options', font=font_settings)
        self.options_frame.grid(row=1, column=0, columnspan=2, sticky='news', padx=xpadding, pady=ypadding)
        
        # convert2cg1 drop down menu
        styles = [True, False]
        self.txtfile = ttk.Combobox(self.options_frame, values=styles, font=font_settings, width=int(maxwidth/10))
        self.txtfile.current(styles.index(txtfile))
        self.txtfile.grid(column=0, row=1)
        self.txtfile_label = tk.Label(self.options_frame, text='txtfile', font=font_settings)
        self.txtfile_label.grid(column=0, row=0)
        
        # N0 entry
        self.N0 = tk.Entry(self.options_frame, font=font_settings, width=int(maxwidth/5))
        self.N0.insert(0, N0)
        self.N0.grid(column=1, row=1)
        self.N0_label = tk.Label(self.options_frame, text='N0 (int)', font=font_settings)
        self.N0_label.grid(column=1, row=0)
        
        # fav entry
        self.fav = tk.Entry(self.options_frame, font=font_settings, width=int(maxwidth/5))
        self.fav.insert(0, fav)
        self.fav.grid(column=2, row=1)
        self.fav_label = tk.Label(self.options_frame, text='fav (float or int)', font=font_settings)
        self.fav_label.grid(column=2, row=0)
        
        # Add padding to all frames in self.options_frame
        for widget in self.options_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))

        #------------#
        # Run button #
        #------------#
        self.run = tk.Button(self.frame, text='Run LUNAR/cluster_analysis.py', font=font_settings, command=self.run_LUNAR)
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
            txt = os.path.join(self.pwd, 'src/GUI_quick_help_pages/cluster_analysis.txt')
            logged = []
            with open(txt, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    # Strip comment's and split by whitespace
                    line = line.split('#')[0]
                    line = line.rstrip()
                    logged.append(line)
        except: # except something failed
            logged.append('FAILED to read LUNAR/src/GUI_quick_help_pages/cluster_analysis.txt document.')
            logged.append('Most likely cause is the cluster_analysis.txt file was renamed in the')
            logged.append('LUNAR/src/GUI_quick_help_pages/ directory or directory names were changed.')
        self.popup(logged, title='Quick help')
        return   
    
    # Function to get filepath for topofile
    def topofile_path(self):
        ftypes = (('data files', '*.data *.data.gz'), ('all files', '*.*'))
        path = filedialog.askopenfilename(initialdir=self.filepath, title='Open topofile?', filetypes=ftypes)
        if path:
            self.filepath = os.path.dirname(os.path.abspath(path))
            path = os.path.relpath(path)
            self.topofile.delete(0, tk.END); self.topofile.insert(0, path);
        return 
    
    # Closing command    
    def closing(self):
        print('Terminating cluster_analysis GUI')
        self.root.destroy()
        return
    
    # Function to run imported LUNAR functions
    def run_LUNAR(self):       
        # Set up log
        log = io_functions.LUNAR_logger()
        valid_inputs = True
        
        # Get information from GUI
        boolean = {'False':False, 'True':True}
        topofile = self.topofile.get()
        try: N0 = int(self.N0.get())
        except:
            log.GUI_error(f'ERROR N0 value of {self.N0.get()} is not an integer value.')
            valid_inputs = False
        try: fav = float(self.fav.get())
        except:
            log.GUI_error(f'ERROR fav value of {self.fav.get()} is not a float value.')
            valid_inputs = False
        txtfile = boolean[self.txtfile.get()]
        
        # Run LUNAR/cluster
        if valid_inputs:
            #------------------#
            # Array processing #
            #------------------#
            if not os.path.isfile(str(topofile)):
                log = io_functions.LUNAR_logger()
                log.configure(level='production', print2console=True, write2log=True)
                files = glob.glob(topofile); array_time = time.time()
                if files:
                    for n, file in enumerate(files, 1):
                        log.clear_all()
                        log.out('\n\nUsing array input option:')
                        log.out(' - topofile     : {}'.format(topofile))
                        log.out(' - matched file : {}'.format(file))
                        log.out(' - elapsed time : {:.2f} (seconds)'.format(time.time() - array_time))
                        log.out(' - progress     : {} of {} ({:.2f}%)'.format(n, len(files), 100*(n/len(files))))
                        try: 
                            inputs = (file, N0, txtfile, fav, log)
                            t1=threading.Thread(target=clusters.analysis, args=inputs)
                            t1.start()
                            t1.join()
                        except Exception:
                            log.GUI_error(traceback.format_exc())
            #------------------#
            # Single usage run #
            #------------------#
            else:
                try: 
                    inputs = (topofile, N0, txtfile, fav, log)
                    t1=threading.Thread(target=clusters.analysis, args=inputs)
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
        N0 = int(self.N0.get())
        fav = float(self.fav.get())
        txtfile = boolean[self.txtfile.get()]
        
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
                if line.startswith('N0') and inputsflag:
                    line = psm.parse_and_modify(line, N0, stringflag=False, splitchar='=')
                if line.startswith('fav') and inputsflag:
                    line = psm.parse_and_modify(line, fav, stringflag=False, splitchar='=')
                if line.startswith('txtfile') and inputsflag:
                    line = psm.parse_and_modify(line, txtfile, stringflag=False, splitchar='=')
                if line.startswith('if __name__ == "__main__":'): inputsflag = False
                f.write(line)
        return