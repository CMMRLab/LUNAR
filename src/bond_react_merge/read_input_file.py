# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 20th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import collections



# Open and read files
def file_inputs(name, log):
    files = {} # { file-tag : path/filename }
    
    # Supported file tags
    supported_filetags = ['pre', 'post', 'data']
    def check_filetags(wholeline, supported_filetags):
        return_boolean = False
        for i in supported_filetags:
            if i in wholeline:
                return_boolean = True
        return return_boolean
    
    # Parser for reading in path and parent_directory
    def parse_in_variables(line):
        string = ''; flag = False; count = 0;
        bounding_characters = ['"', "'"] # " or ' bounding chars
        for i in line:
            if i in bounding_characters: flag = True; count += 1; continue;
            if count >= 2: flag = False
            if flag: string += i
        return string.rstrip()
    
    # open file to read inouts
    with open(name, 'r') as f:
        # Inialize flags
        path_flag = False
        file_flag = False
        dir_flag = False
        class_flag = False
        
        # set intial default values
        path = ''              # Initialize as empty and update if found
        parent_directory = ''  # Initialize as empty and update if found
        ff_class = False       # Initialize as False and update if found
        tmpfiles = {}          # { file-tag : path/filename }
        filetags = []          # lst to check for duplicate filetags
        
        # loop through each line in file
        for wholeline in f:
            
            # Strip comment's
            wholeline = wholeline.split('#')[0]
            wholeline = wholeline.rstrip()

            # strip and split line by white space                
            line = wholeline.strip()
            line = line.split()

            # Set flags
            if line == [] or '#' in wholeline:
                path_flag = False
                file_flag = False
                dir_flag = False
                class_flag = False
            elif 'path' in wholeline:
                path_flag = True
                file_flag = False
                dir_flag = False
                class_flag = False
            elif 'path' in wholeline:
                path_flag = True
                file_flag = False
                dir_flag = False
                class_flag = False
            elif 'parent_directory' in wholeline:
                dir_flag = True
                path_flag = False
                file_flag = False
                class_flag = False
            elif 'ff_class' in wholeline:
                file_flag = False
                path_flag = False
                dir_flag = False
                class_flag = True
            elif check_filetags(wholeline, supported_filetags) and len(line) == 2:
                file_flag = True
                path_flag = False
                dir_flag = False
                class_flag = False
                
            # Find path
            if path_flag:
                # Parse wholeline to get info
                tmpline = parse_in_variables(wholeline)
                
                # Check if tmpline begins with '/' or '\' character if so remove
                if len(tmpline) > 0:
                    if tmpline[0] == '\\': tmpline = tmpline.replace('\\', '', 1) 
                    if tmpline[0] == '/': tmpline = tmpline.replace('/', '', 1) 
                
                # Check if tmpline ends with '/' character if not add to tmpline
                if len(tmpline) > 0:
                    if not tmpline[-1] == '/': tmpline = '{}/'.format(tmpline)
                
                # Update path string
                path = tmpline
                
            # Find dir
            if dir_flag:
                # Parse wholeline to get info
                tmpline = parse_in_variables(wholeline)
                
                # Check if tmpline begins with '/' or '\' character if so remove
                if len(tmpline) > 0:
                    if tmpline[0] == '\\': tmpline = tmpline.replace('\\', '', 1) 
                    if tmpline[0] == '/': tmpline = tmpline.replace('/', '', 1) 
                
                # Check if tmpline ends with '/' character if not add to tmpline
                if len(tmpline) > 0:
                    if not tmpline[-1] == '/': tmpline = '{}/'.format(tmpline)
                    
                # Update path string
                parent_directory = tmpline
                
            # Find ff class
            if class_flag: ff_class = int(wholeline[-1])

            
            # Find files
            elif file_flag:
                filetag = line[0]
                filename = line[1]
                tmpfiles[filetag] = filename
                filetags.append(filetag)

        # Check for duplicate filetags
        duplicates = [item for item, count in collections.Counter(filetags).items() if count > 1]
        if duplicates:
            for i in duplicates:
                log.error(f'ERROR filetag: {i} is used more then once in inputfile {name}')
            
        # Create files dictionary with path append to file name
        for i in tmpfiles:
            files[i] = '{}{}'.format(path, tmpfiles[i])
            
        # Print the read in files to the user
        log.out(f'\n\nThe following files and their corresponding filetags have been read in from {name}:')
        for i in files:
            log.out('{:^8} {:^40}'.format(i, files[i]))
        
        # Print the parent_directory
        if parent_directory != '':
            log.out(f'parent_directory was set from {name} as {parent_directory}')
        return files, parent_directory