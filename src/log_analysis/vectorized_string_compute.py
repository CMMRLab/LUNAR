# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
February 13th, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

    *********************************************************
    * Requirements:                                         *
    *   python numpy module:                                *
    *    - pip3 install numpy (if pip manager is installed) *
    *********************************************************
    
https://www.geeksforgeeks.org/eval-in-python/
"""
##############################
# Import Necessary Libraries #
##############################
import numpy as np
import traceback



##########################################################################
# Function to compute vectorized operations on arrays using string math. #
#  Inputs:                                                               #
#     data = {'NAME':[lst], ... }                                        #
#     compute_string = "${A} + ${B}/2 + 100", where A and B are keys in  #
#                       the data dictionary.                             #
#  Outputs:                                                              #
#    computed_list = "${A} + ${B}/2 + 100", using vectorization methods  #
#                    from numpy and Python eval() function.              #
##########################################################################
def thermo_data(data, compute_string):
    #--------------------------------------#
    # d will be used in list(eval(string)) #
    #--------------------------------------#
    d = {i:np.array(data[i]) for i in data if isinstance(data[i], tuple) or isinstance(data[i], list)}
    if d: nzeros = max([len(d[i]) for i in d])
    else: nzeros = 1
    message = ''
    
    #------------------------------------------------------------------#
    # The compute_string needs to be transformed, for example:         #
    #  in-string:  "  ${A} + ${B}/2   + 100"                           #
    #  out-string: "d['A'] + d['B']/2 + 100"                           #
    # where A and B are columns in the logfile, stored as keys in the  #
    # dictionary d = {'COLUMN-NAME':np.array()}. Using the Python      #
    # eval() function we can then perform vectorized computes of the   #
    # numpy arrays.                                                    #
    #------------------------------------------------------------------#
    cstring = compute_string
    cstring = cstring.replace("${", "d['")
    cstring = cstring.replace("}", "']")
    
    #---------------------------------------------------------------#
    # Check that each variable in the compute_string exists in the  #
    # dictionary "d". Warn if not in the dictionary.                #
    #---------------------------------------------------------------#
    # Find variables from compute_string
    variables = set(); tmp = ''; dollar_sign = False; bracket0 = False; bracket1 = False;
    for i in compute_string:
        if i == '$': dollar_sign = True
        if i == '{': bracket0 = True
        if i == '}': bracket1 = True
        if dollar_sign and bracket0:
            tmp += str(i)
        if bracket1:
            dollar_sign = False; bracket0 = False; bracket1 = False;
            tmp = tmp.replace('{', '')
            tmp = tmp.replace('}', '')
            variables.add(tmp); tmp = '';
            
    # Check that each variable is in the "d" dictionary
    variable_check = True
    for variable in variables:
        if variable not in d:
            variable_check = False
            message += 'ERROR {} is not in the compute dictionary. Perhaps something was misspelled.\n'.format(variable)     
        
    #---------------------------------------------------------------------#
    # Try making compute, if it fails log a different message and create  #
    # an array of zeros to output.                                        #
    #---------------------------------------------------------------------#
    if variable_check:
        # Setup the globals namespace to limit scope of what eval() can do
        allowed_builtins = ['min','max','sum','abs','len','map','range','reversed']
        copied_builtins = globals()['__builtins__'].copy()
        globals_dict = {'d':d}
        globals_dict['__builtins__'] = {key:copied_builtins[key] for key in allowed_builtins}
        globals_dict['np'] = globals()['np']
    
        # Try computing
        try:
            new_data = list(eval(cstring, globals_dict))
            message += f'"{compute_string}" was succesfully evaluted.\n'
        except Exception:
            new_data = nzeros*[0]
            message += f'ERROR could not evalute "{compute_string}" compute.\n'
            message += '\n{}\n'.format(str(traceback.format_exc()))
    else: new_data = nzeros*[0]
    return new_data, message


if __name__ == "__main__":
    data = {'test':'loop '}
    compute_string = "__import__('os').getlogin()"
    #compute_string = "[d['test'] +let for let in 'abc']"
    x = thermo_data(data, compute_string)
    print(x[0])
    print(x[1])