# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
March 28th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import sys
import os


##################################################################################
# This class is meant to handle LUNAR's logging. If you are going to add onto    #
# LUNAR, it is recommended to make use of this logging class to stay consistent  #
# with the how LUNAR, log's results. In general follow these rules:              #
#   Generate an instance of the class as "log"                                   #
#      log = LUNAR_log()                                                         #
#   How to log using methods                                                     #
#     log.out('string'), where 'string' can be a f-string. Use instead of print  #
#     log.warn('string'), where 'string' can be a f-string and the should start  #
#                         start with WARNING (i.e. 'WARNING inconsistent ...')   #
#     log.error('string'), where 'string' can be a f-string and should start     #
#                          with ERROR (i.e. 'ERROR filename does not exist ...').#
#                          the error method will also call sys.exit() and will   #
#                          kill the program so only use for critically incorrect #
#                          checks.                                               #
#                                                                                #
# This logging structure was generated after some decent thought on how things   #
# should be logged and it was decided to not use the standard python logging     #
# module. The main purpose of this logging class is to capture the console       #
# outputs, when using the GUI run modes, such that the console outputs can be    #
# written to a pop-up window. Every main function in LUNAR's modules will have   #
# the default log class set as def main(topofile, ..., log=LUNAR_log()): where   #
# each module will have the log class encapsoluted into the main function. The   #
# the log class could also be instantiated, before calling the main function as  #
#   log = LUNAR_log()                                                            #
#   main(topofile, ..., log=log)                                                 #
#   print(''.join(log.logged))                                                   #
# where all logged strings can be accessed after the main function has completed #
# or exited.                                                                     #
##################################################################################
class LUNAR_logger:
    def __init__(self, level='production', print2console=True, write2log=True):
        self.logged = []   # list of strings that have been logged (outputs, warnings, and errors)
        self.warnings = [] # list of strings from warnings
        self.errors = []   # list of strings from errors
        self.debugs = []   # list of strings from debugs
        self.level = level
        self.write2log = write2log
        self.print2console = print2console
        
        # This scripts path and name
        script_path = os.path.dirname(os.path.realpath(__file__))
        #script_name = os.path.basename(__file__)
        self.lunar_path = os.path.normpath(os.path.join(script_path, '../'))
        
    # Use this method instead of "print()" command
    def out(self, text):
        if self.print2console: print(text)
        self.logged.append(text)
        
    # Use this method to print and log warnings
    def warn(self, text):
        if self.print2console: print(text)
        self.logged.append(text)
        self.warnings.append(text)
        
    # Use this method to print and log errors
    def error(self, text):
        if self.print2console: print(text)
        self.logged.append(text)
        self.errors.append(text)
        self.out('Changing back to: "{}" directory before exiting'.format(self.lunar_path))
        os.chdir(self.lunar_path)
        sys.exit()
        
    # use this method for a GUI error (not exiting)
    def GUI_error(self, text):
        if self.print2console: print(text)
        self.logged.append(text)
        self.errors.append(text)
        
    # use this method for debugging code
    def debug(self, text):
        if self.level == 'debug':
            debug_text = '{}  DEBUG'.format(text)
            if self.print2console: print(debug_text)
            self.debugs.append(text)
            self.logged.append(debug_text)
        
    # Use this method to show number of warnings and errors
    def out_warnings_and_errors(self):
        self.out('\n\nNumber of WARNING(s) and ERROR(s)')
        self.out('  {} WARNING(s)'.format(len(self.warnings)))
        self.out('  {} ERROR(s)'.format(len(self.errors)))
        
    # use this method to clear all entries in logged, warnings, and errors
    def clear_all(self):
        self.logged = []; self.warnings = []; self.errors = [];  
        
    # use this method to write all logged values (let parent_directory take care of location)
    def write_logged(self, filename):
        if self.write2log:
            with open(filename, 'w') as f:
                for text in self.logged:
                    f.write('{}\n'.format(text))

    # use this method to configure logger from defaults
    def configure(self, level='', print2console='', write2log=''):
        # Update level
        if level in ['production', 'debug'] and level != self.level:
            #print('Swithing level of logger to {}.'.format(level))
            self.level = level
            
        # update print2console
        if print2console in [True, False] and print2console != self.print2console:
            #print('Swithing print2console of logger to {}.'.format(print2console))
            self.print2console = print2console
            
        # update write2log
        if write2log in [True, False] and write2log != self.write2log:
            #print('Swithing write2log of logger to {}.'.format(write2log))
            self.write2log = write2log


############################
# Function to get basename #
############################
def get_basename(topofile, newfile='', character=':', pflag=True, log=LUNAR_logger()):
    root = os.path.basename(topofile)
    if root.endswith('data.gz'): 
        root = root.replace('.gz', '') # remove .gz suffix
    root = root[:root.rfind('.')]
    character_check = newfile.strip() # avoid whitespace errors
    if character_check.startswith(character): # prefix
        basename = '{}{}'.format(root, newfile[1:])
    elif character_check.endswith(character): # suffix
        basename = '{}{}'.format(newfile[:-1], root)
    elif newfile in ['', ' ', '   ']: # single, double, or triple spacing
        if pflag:
            log.warn('WARNING newfile is an empty string. This may result in overwritting original filenames!')
            log.out(f'  Please be aware that any file with basename {root}')
            log.out('  maybe overwritten.')
        basename = '{}'.format(root)
    else: # reset completely
        basename = '{}'.format(newfile)
    return basename


###################################################
# Function to build directories based on topofile #
###################################################
def get_dir_from_topofile(topofile, parent_directory):
    if '/' in parent_directory:
        split = parent_directory.split('/')
    elif '\\' in parent_directory:
        split = parent_directory.split('\\')
    else: split = []
    if len(split) > 0:
        if split[0] == 'topofile':
            added = '/'.join(split[1:])
            path = os.path.dirname(os.path.abspath(topofile))
            path = os.path.join(path, added)
        else: path = os.path.dirname(os.path.abspath(topofile))
    else: path = os.path.dirname(os.path.abspath(topofile))
    return path


###############################################################
# When python loads a path it may set the directory delimeter #
# as a '\' however this is not correct for setting in a file  #
# so this function converts  any '\' to '/' characters        #
###############################################################
def path_to_string(path):
    string = str(path)
    string = string.replace('\\', '/')
    return string


#############################################################
# Function to check if topofile, logfile, or xrdfile or ... #
# already exists in a given parent_directory to avoid re-   #
# processing in a "glob" array processing run.              #
#############################################################
def check_outfile_existance(filename, newfile, parent_directory, filetype='topofile', check_per_file=True):
    # Default is to check each file's existance, prior to processing. If users don't like this
    # behavior all they have to do is come to this function and set the check_per_file=False
    # which will then change this behavior across all LUNAR array processing runs.
    if check_per_file:
        # Get file basename
        if filetype != 'topofile': filename = filename.replace(filetype, 'topofile')
        basename = get_basename(filename, newfile=newfile, character=':', pflag=True)
        
        # Find present working directory
        pwd = os.getcwd()
        
        # Find/create paths to store code results
        path = os.path.join(pwd, parent_directory)
        
        # If parent_directory == 'topofile' use topofile path as path
        if filetype in parent_directory:
            if filetype != 'topofile': parent_directory = parent_directory.replace(filetype, 'topofile')
            path = get_dir_from_topofile(filename, parent_directory)
        
        # Set full file path to the output '*.log.lunar' file
        full_path = os.path.join(path, '{}.log.lunar'.format(basename))
        try: file_exists = os.path.isfile(full_path)
        except: file_exists = False
    else: file_exists = False
    return file_exists