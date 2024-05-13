# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
February 16, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

    ***************************************************
    * Requirements:                                   *
    *   python 3.7+                                   *
    *                                                 *
    ***************************************************
"""

############################################################
# Helper functions to use outside of this file for getting #
# sections from thermo file (appends zeros when needed)    #
############################################################
# Function to check if value is an int
def check_int(value):
    try: 
        int(value)
        return True
    except: 
        return False
    
# Function to get most frequent element in list
def most_frequent(List):
    return max(set(List), key = List.count)

# Function to get int valued sectionIDs from strings like:
#    'all' or  '1,2,3'  or  '1-4'  or  '1,3-5'
def get_sections(log_sections, sections):
    sectionIDs = []
    if 'all' in sections:
        sectionIDs = log_sections
    else:
        IDs = sections.split(',')
        for i in IDs:
            if '-' in i:
                span = i.split('-')
                if len(span) == 2:
                    if check_int(span[0]) and check_int(span[1]):
                        span = [int(j) for j in span]
                        lo = min(span); hi = max(span);
                        for j in range(hi - lo + 1):
                            if j+lo in log_sections:
                                sectionIDs.append(j+lo)
                    else: print(f'  WARNING span {i} value(s) are not an int values. skipping.')
                else: print(f'  WARNING section span was not of length two {i}')
            else:
                if check_int(i):
                    if int(i) in log_sections:
                        sectionIDs.append(int(i))
                else: print(f'  WARNING ID {i} is not an int. skipping.')
    return sorted(set(sectionIDs))

# Function to get all sections that user desires from log file
def get_data(log, sections, pflag=True):
    # Find sectionIDs based on user specified string
    sectionIDs = get_sections(log.sections, sections)
    if pflag: print(f'  Loading {" ".join([str(i) for i in sectionIDs])} of {" ".join([str(i) for i in log.sections])} sections')
    
    # Find all column names in loaded sections and start joining data
    columns = {column for ID in log.sections for column in log.data[ID]}
    data = {} # {column-name:[lst of data]}
    for ID in sectionIDs:
        used = {column:False for column in columns}; ndatas = [];
        for column in log.data[ID]:
            used[column] = True
            ndatas.append(len(log.data[ID][column]))
            if column in data:
                data[column].extend(log.data[ID][column])
            else: data[column] = log.data[ID][column]

        # Check for any unused columns and make zeros to append
        for column in used:
            if not used[column]:
                zeros = [0]*int(most_frequent(ndatas))
                if column in data:
                    data[column].extend(zeros)
                else: data[column] = zeros
    return data


####################################################
# Functions for parsing and converting thermo data #
####################################################
# Function to count number of matching keywords per line
def count_keywords(line, keywords):
    count = 0
    for i in line:
        if i in keywords: count += 1
    return count

# Function to convert all strings in line to floats of ints
def line2digits(line, nline):
    digits = []
    for i in line:
        # 1st try converting to float, 2nd check if float is really and not int, and 3rd log digit as int or float
        try: digit = float(i)
        except: raise Exception(f'  ERROR not all thermo data could not be converted to floats or ints at line {nline}\n{line}')
        if digit.is_integer(): digit = int(digit)
        digits.append(digit)
    return digits


#####################################################################
# Class to open and parse LAMMPS log file into data structures.     #
#   Inputs:                                                         #
#      logfile  = name of log file to read                          #
#      keywords = list of keywords for column data, default is      #
#                 'Step' only. Add keywords when necessary.         #
#   Outputs:                                                        #
#      A class with the following attributes:                       #
#         CLASS.data = { sectionID : {column1 : [column1], ... } }  #
#         CLASS.sections  = [ lst of section IDs ]                  #
#         CLASS.sectionID = count of sections in log file           #
#         CLASS.nlines = count of lines in file                     # 
#         CLASS.nentries = count of entries in file                 #
#####################################################################
class file:
    def __init__(self, logfile, keywords=['Step'], pflag=True):
        self.data = {} # { sectionID : {column1 : [column1], ... } }
        self.sections = [] # [ lst of section IDs ]
        self.sectionID = 0 # tally of section count
        self.nlines = 0 # count of lines in file
        self.nentries = 0 # count of entries in file
        
        # Flags to find things
        data_flag = False
        
        # Open and parse log file
        with open(logfile, 'r') as f:
            for string in f:
                # Strip and split line
                line = string.strip(); line = line.split(); self.nlines += 1;
                
                # Find begining thermo data section if count of keywords is > 0
                if 'WARNING' in string or string == '': continue # skip warnings or empty lines
                elif count_keywords(line, keywords) > 0 and '#' not in string and 'thermo_style' not in string: # section of logfile with thermo data
                    data_flag = True
                    self.sectionID += 1
                    self.sections.append(self.sectionID)
                    self.data[self.sectionID] = {i:[] for i in line} # {column-name:[]}
                    indexes = {n:i for n, i in enumerate(line)} # {index:column-name}
                    continue
                # Break thermo data section if 'Loop time of' in string
                elif 'Loop time of' in string: data_flag = False
                elif 'Loop' in string and 'time' in string and 'of' in string: data_flag = False
                elif 'ERROR' in line: data_flag = False
                    
                # Add data to self.data if data_flag
                if data_flag:
                    digits = line2digits(line, self.nlines) # convert strings to numbers (ints and floats)
                    if len(digits) == len(indexes): # skip over incomplete lines
                        for n, i in enumerate(digits):
                            self.data[self.sectionID][indexes[n]].append(i)
                            self.nentries += 1
                    else: 
                        print(f'  WARNING skipping over line {self.nlines} in log file due to incomplete data series. Line {self.nlines}: ')
                        print('     {}'.format( ' '.join(line) ))
        
        # Print to user to outcomes
        if pflag:
            print(f'  read sections: {" ".join([str(i) for i in self.sections])}')
            print(f'  with {self.nentries} log entries')