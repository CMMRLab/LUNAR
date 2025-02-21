# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
February 18, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import difflib
import glob
import os


# Function to convert string to float or int and will default to string if needed
def string2digit(string):
    digit = string
    try: 
        digit = float(string)
        if digit.is_integer():
            digit = int(digit)
    except: pass
    return digit


# Function to get a value bound between two characters in a string. Example:
#    value = get_value_between_characters("example[EXTRACT-ME]example")
#    value = EXTRACT-ME
def get_value_between_characters(string, start_char="[", end_char="]"):
    start_index = string.find(start_char)
    if start_index == -1:
        return None

    end_index = string.find(end_char, start_index + 1)
    if end_index == -1:
        return None
    return string2digit(string[start_index + 1:end_index])

    
# Function to get wildcards between a glob_string and a string glob matched with. Inspired from:
# https://stackoverflow.com/questions/17904097/python-difference-between-two-strings
def get_glob_wildcards(glob_string, matched_string):
    # Force both strings to follow the same format rules
    glob_string = os.path.normpath(glob_string)
    matched_string = os.path.normpath(matched_string)
    
    # Find all strings that diff sets as adds "+"
    wildcards = []; add = False; tmp = ''
    for i, s in enumerate(difflib.ndiff(glob_string, matched_string)):
        if s[0] == ' ':
            add = False
        if s[0] == '-':
            if tmp != '': wildcards.append(tmp)
            add = False; tmp = ''
        if s[0] == '+':
            add = True 
        if add:
            tmp += s[-1]
    if tmp != '': wildcards.append(tmp)
            
    # Check to see if any wildcards show up in glob_string, if they do, remove them
    # from wildcards list. Additionally try converting from a string to float or int
    wildcards = [string2digit(i) for i in wildcards if i not in glob_string]
    return wildcards


# Function to find wildcards that are "most numeric" and convert
# them to a float or int value. Useful for atom_typing log analysis
def get_numeric_wildcard(wildcards):
    numeric = None # default return incase there are no wildcards or wildcard is all alpha-characters
    if wildcards: 
        # Find the %number of characters that are digits in each wildcard
        digitized = {} # {index:(wildcard, %digits), ... nwildcards}
        for n, wildcard in enumerate(wildcards):
            wildcard = str(wildcard); ndigits = 0
            for i in wildcard:
                try: 
                    float(i)
                    ndigits += 1
                except: pass
            percent_digits = 100*(ndigits/len(wildcard))
            digitized[n] = (wildcard, percent_digits)
            
        # Select the wild card with the largest percent digits as the
        # numeric value to return to use as the independant variable
        maximized_digits = max(digitized.items(), key=lambda x: x[1][1] ) # [0=keys;1=values][1=index in value tuple]
        maximized_wildcard = maximized_digits[1][0]; digits = ''
        for i in maximized_wildcard:
            try: 
                float(i)
                digits += i
            except: pass
        if digits:
            numeric = float(digits)
            if numeric.is_integer():
                numeric = int(numeric)
    return numeric


# Basic test of the log file reader             
if __name__ == "__main__":    
    print('\n')
    print('+------------+')
    print('|   Test 1   |')
    print('+------------+')
    topofile = '../../../**EXAMPLES/array_processing/*.[mp][od][lb]*'
    topofiles = glob.glob(topofile)
    for file in topofiles:
        josh_wildcards = get_glob_wildcards(topofile, file)
        numeric = get_numeric_wildcard(josh_wildcards)
        print(topofile)
        print(file)
        print('josh_wildcards = ', josh_wildcards)
        print('  numeric = ', numeric)
        print('  type(numeric) = ', type(numeric))
        print()
    
    print('\n')
    print('+------------+')
    print('|   Test 2   |')
    print('+------------+')
    topofile = '../../../EXAMPLES/array_processing/?????.*'
    topofiles = glob.glob(topofile)
    for file in topofiles:
        josh_wildcards = get_glob_wildcards(topofile, file)
        numeric = get_numeric_wildcard(josh_wildcards)
        print(topofile)
        print(file)
        print('josh_wildcards = ', josh_wildcards)
        print('  numeric = ', numeric)
        print('  type(numeric) = ', type(numeric))
        print()
        
    print('\n')
    print('+------------+')
    print('|   Test 3   |')
    print('+------------+')
    topofile = '../../../EXAMPLES/array_processing/poly_tracking_replicate_1_time_*.data'
    topofiles = glob.glob(topofile)
    for file in topofiles:
        josh_wildcards = get_glob_wildcards(topofile, file)
        numeric = get_numeric_wildcard(josh_wildcards)
        print(topofile)
        print(file)
        print('josh_wildcards = ', josh_wildcards)
        print('  numeric = ', numeric)
        print('  type(numeric) = ', type(numeric))
        print()
        
    print('\n')
    print('+------------+')
    print('|   Test 4   |')
    print('+------------+')
    topofile = '../../../EXAMPLES/array_processing/poly_tracking_replicate_1_time_*ps.data'
    topofiles = glob.glob(topofile)
    for file in topofiles:
        josh_wildcards = get_glob_wildcards(topofile, file)
        numeric = get_numeric_wildcard(josh_wildcards)
        print(topofile)
        print(file)
        print('josh_wildcards = ', josh_wildcards)
        print('  numeric = ', numeric)
        print('  type(numeric) = ', type(numeric))
        print()
     
    print('Testing how to get index from string')
    s = 'wildcards[0]'
    print(' -', s)