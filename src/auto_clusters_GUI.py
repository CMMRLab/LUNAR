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
from tkinter.scrolledtext import ScrolledText
from tkinter import filedialog
from tkinter import Toplevel
from tkinter import ttk
import tkinter as tk
import traceback
import math
import sys
import os


# setting path to files called by other codes in src
sys.path.append('../')
from auto_cluster_analysis import main


###################################
# LUNAR/auto_cluster_analysis GUI #
###################################
class GUI:
    def __init__(self, files_directory, N0, fav, txtfile, newfile, GUI_zoom):
        
        # Find present working directory
        self.pwd = os.getcwd()
        self.filepath = self.pwd
        self.filename = os.path.join(self.pwd, 'auto_cluster_analysis.py')

        # Initialize window
        self.root = tk.Tk()
        self.root.title('LUNAR/auto_cluster_analysis.py GUI v1.0')
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
        self.maxwidth = 75
        
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
        

        # files_directory entry
        self.files_directory = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
        self.files_directory.insert(0, files_directory)
        self.files_directory.grid(column=1, row=4)
        self.dir_button = tk.Button(self.inputs_frame, text='files_directory', font=font_settings, command=self.directory_path)
        self.dir_button.grid(column=0, row=4)
        
        # newfile entry
        self.newfile = tk.Entry(self.inputs_frame, width=maxwidth, font=font_settings)
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
        self.options_frame.grid(row=1, column=0, sticky='news', padx=xpadding, pady=ypadding)
        
        # convert2cg1 drop down menu
        styles = [True, False]
        self.txtfile = ttk.Combobox(self.options_frame, values=styles, width=int(maxwidth/7.5), font=font_settings)
        self.txtfile.current(styles.index(txtfile))
        self.txtfile.grid(column=0, row=1)
        self.txtfile_label = tk.Label(self.options_frame, text='txtfile', font=font_settings)
        self.txtfile_label.grid(column=0, row=0)
        
        # N0 entry
        self.N0 = tk.Entry(self.options_frame, width=int(maxwidth/3.75), font=font_settings)
        self.N0.insert(0, N0)
        self.N0.grid(column=1, row=1)
        self.N0_label = tk.Label(self.options_frame, text='N0 (int)', font=font_settings)
        self.N0_label.grid(column=1, row=0)
        
        # fav entry
        self.fav = tk.Entry(self.options_frame, width=int(maxwidth/3.75), font=font_settings)
        self.fav.insert(0, fav)
        self.fav.grid(column=2, row=1)
        self.fav_label = tk.Label(self.options_frame, text='fav (float or int)', font=font_settings)
        self.fav_label.grid(column=2, row=0)
        
        # Add padding to all frames in self.inputs_frame
        for widget in self.options_frame.winfo_children():
            widget.grid_configure(padx=xpadding, pady=int(ypadding/2))

        #------------#
        # Run button #
        #------------#
        self.run = tk.Button(self.frame, text='Run LUNAR/auto_cluster_analysis.py', font=font_settings, command=self.run_LUNAR)
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
    # Function to get directory
    def directory_path(self):
        path =filedialog.askdirectory(initialdir=self.filepath)
        if path:
            self.filepath = os.path.dirname(os.path.abspath(path))
            path = os.path.relpath(path)
            self.files_directory.delete(0, tk.END); self.files_directory.insert(0, path);
        return
    
    # Closing command    
    def closing(self):
        print('Terminating auto_cluster_analysis GUI'); self.root.destroy();
        return
    
    # Function to run imported LUNAR functions
    def run_LUNAR(self):  
        # Set up log
        log = io_functions.LUNAR_logger()
        valid_inputs = True
        
        # Get information from GUI
        boolean = {'False':False, 'True':True}
        files_directory = self.files_directory.get()
        newfile = self.newfile.get()
        try: N0 = int(self.N0.get())
        except:
            log.GUI_error(f'ERROR N0 value of {self.N0.get()} is not an integer value')
            valid_inputs = False
        try: fav = float(self.fav.get())
        except:
            log.GUI_error(f'ERROR fav value of {self.fav.get()} is not a float value')
            valid_inputs = False
        txtfile = boolean[self.txtfile.get()]
        
        # Run LUNAR/auto_cluster_analysis
        if valid_inputs:
            try: main(files_directory, N0, fav, txtfile, newfile, log=log)
            except Exception:
                log.GUI_error(traceback.format_exc())
        self.popup(log.logged)
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
        files_directory = io_functions.path_to_string(self.files_directory.get())
        newfile = self.newfile.get()
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
                if line.startswith('files_directory') and inputsflag:
                    line = psm.parse_and_modify(line, files_directory, stringflag=True, splitchar='=')
                if line.startswith('newfile') and inputsflag:
                    line = psm.parse_and_modify(line, newfile, stringflag=True, splitchar='=')
                if line.startswith('N0') and inputsflag:
                    line = psm.parse_and_modify(line, N0, stringflag=False, splitchar='=')
                if line.startswith('fav') and inputsflag:
                    line = psm.parse_and_modify(line, fav, stringflag=False, splitchar='=')
                if line.startswith('txtfile') and inputsflag:
                    line = psm.parse_and_modify(line, txtfile, stringflag=False, splitchar='=')
                if line.startswith('if __name__ == "__main__":'): inputsflag = False
                f.write(line)
        return