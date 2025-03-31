# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
June 4th, 2024
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
    
# Function to check if value is a float
def check_float(value):
    try: 
        float(value)
        return True
    except: 
        return False
    
# Function to get most frequent element in list
def most_frequent(List):
    return max(set(List), key = List.count)

# Function to get int valued sectionIDs from strings like:
#    'all' or  '1,2,3'  or  '1-4'  or  '1,3-5'
def get_sections(log_sections, sections):
    # Force sections to be a string for this, so users can supply ints if they want
    sections = str(sections); sectionIDs = []
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


####################################################
# Functions for parsing and converting thermo data #
####################################################
# Function to count number of matching keywords per line
def count_keywords(line, keywords):
    count = 0
    for keyword in keywords:
        count += line.count(keyword)
    return count

# Function to convert string to float or int and will default to string if needed
def string2digit(string):
    digitized = True; digit = string
    try: 
        digit = float(string)
        if digit.is_integer():
            digit = int(digit)
    except: digitized = False
    return digitized, digit


######################################################################
# Class to open and parse LAMMPS log file into data structures.      #
#   Inputs:                                                          #
#      logfile  = name of log file to read                           #
#      keywords = list of keywords for column data, default is       #
#                 'Step' only. Add keywords when necessary.          #
#   Outputs:                                                         #
#      A class with the following attributes:                        #
#         CLASS.data = { sectionID : {column1 : [column1], ... } }   #
#         CLASS.sections  = [ lst of section IDs ]                   #
#         CLASS.sectionID = count of sections in log file            # 
#         CLASS.nlines = count of lines in file                      # 
#         CLASS.nentries = count of entries in file                  #
#                                                                    #
#   Methods:                                                         #
#      .get_data(self, sections, remove_duplicates=True, pflag=True) #
#         Which allows for multiple sections of a logfile to be      #
#         "joined" together based on sections string.                #
######################################################################
class file:
    def __init__(self, logfile, keywords=['Step'], pflag=True):
        self.data = {} # { sectionID : {column1 : [column1], ... } }
        self.sections = [] # [ lst of section IDs ]
        self.sectionID = 0 # tally of section count
        self.nlines = 0 # count of lines in file
        self.nentries = 0 # count of entries in file
        
        # nmax for printing sections with pflag
        self.nmax = 10
        
        # Open and parse log file
        data_flag = False
        with open(logfile, 'r') as f:
            for string in f:
                line = string.strip(); line = line.split(); self.nlines += 1;
                if 'WARNING' in string or len(string) == 0: continue
                elif count_keywords(line, keywords) > 0 and '#' not in string and 'thermo_style' not in string:
                    data_flag = True; self.sectionID += 1
                    self.sections.append(self.sectionID) 
                    self.data[self.sectionID] = {i:[] for i in line}
                    indexes = {n:i for n, i in enumerate(line)}; continue
                elif not line: data_flag = False
                elif 'Loop time of' in string: data_flag = False
                elif 'Loop' in string and 'time' in string and 'of' in string: data_flag = False
                elif 'ERROR' in line: data_flag = False
                if data_flag:
                    if len(line) != len(indexes):
                        print(f'\n  WARNING skipping over line {self.nlines} in log file due to incomplete data series. Line {self.nlines}: ')
                        print('     {}'.format(' '.join(line) ))
                    else:
                        success = True
                        for n, i in enumerate(line):
                            digitized, digit = string2digit(i)
                            if digitized:
                                column = indexes[n]; self.nentries += 1
                                self.data[self.sectionID][column].append(digit)
                            else: success = False; break
                        if not success:
                            print(f'\n  WARNING not all thermo data could not be converted to floats or ints at line {self.nlines}:')
                            print('     {}'.format(' '.join(line) ))
                            for i in range(n):
                                column = indexes[i]; self.nentries -= 1
                                del self.data[self.sectionID][column][-1]
                                

        # Print to user to outcomes
        if pflag:
            if len(self.sections) < self.nmax:
                print(f'  read sections: {" ".join([str(i) for i in self.sections])}')
            print(f'  with {self.nentries} log entries')
            
    # method to get all sections that user desires from log file
    def get_data(self, sections, remove_duplicates=True, pflag=True):
        # Find sectionIDs based on user specified string
        sectionIDs = get_sections(self.sections, sections)
        if pflag and len(self.sections) < self.nmax:
            print(f'  Loading {" ".join([str(i) for i in sectionIDs])} of {" ".join([str(i) for i in self.sections])} sections')
        
        # Find all column names in loaded sections and start joining data
        columns = {column for ID in self.sections for column in self.data[ID]}
        data = {} # {column-name:[lst of data]}
        for ID in sectionIDs:
            used = {column:False for column in columns}; ndatas = [];
            for column in self.data[ID]:
                used[column] = True
                tmp = self.data[ID][column]
                if column in data:
                    if remove_duplicates and data[column]:
                        try: del tmp[0]
                        except: pass
                    data[column].extend(tmp)
                else: data[column] = tmp
                ndatas.append(len(tmp))
    
            # Check for any unused columns and make zeros to append
            for column in used:
                if not used[column]:
                    zeros = [0]*int(most_frequent(ndatas))
                    if column in data:
                        data[column].extend(zeros)
                    else: data[column] = zeros
        return data