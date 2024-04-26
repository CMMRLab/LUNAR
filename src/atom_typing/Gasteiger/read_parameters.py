# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
February 21st, 2022
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Function to check if variable is a float
def check_float(variable):
    try:
        float(variable)
        return_boolean = True
    except:
        return_boolean = False
    return return_boolean


# class to create ojbect
class Parameters:
    pass # .a .b. c 


# Open and read the file into a dictionary
def read(filename):
    parameters = {} # { element : {hybridization tag : Parameters Object} }

    # Open and read file
    with open(filename, 'r') as f:
        
        # Initializing flags    
        parms_flag = False
    
        # Loop through each line
        for line in f:
            
            # Strip comment's
            line = line.split('!')[0]
            line = line.rstrip()
            
            # Strip line
            line = line.strip()
            line = line.split()
            
            # Setting flags if '#' character is in string of 1st element in list
            # as False to break each section flag from previous iteration
            if len(line) > 0 and line[0][0] == '#':
                parms_flag = False
                
            # parameters info
            if '#parameters' in line:
                parms_flag = True
                
            # Read parameters
            if parms_flag:
                # Check if len(line) >= 5 and if line[2, 3, 4] is a float
                if len(line) >= 5 and check_float(line[2]) and check_float(line[3]) and check_float(line[4]):
                    #print(line)
                    element = line[0]
                    hybridization = line[1]
                    
                    # Set parameters object
                    P = Parameters()
                    P.a = float(line[2])
                    P.b = float(line[3])
                    P.c = float(line[4])
                    
                    # Add object if element already in parameters
                    if element in parameters:
                        parameters[element][hybridization] = P
                    else:
                        parameters[element] = {hybridization : P}

    return parameters

                    
                    
#############################    
### Testing reading parms ###    
#############################
if __name__ == '__main__':
    
    # frc FF file to read
    file = 'Gasteiger_parameters.txt'    
    parms = read(file)
    
    # Try accessing some parms
    
    # Sp1 Carbon
    parm = parms['C']['Sp1']
    print('\n\n{:^24}'.format('Sp1 Carbon'))
    print('{:^8} {:^8} {:^8}'.format('a', 'b', 'c'))
    print('{:^8} {:^8} {:^8}'.format(parm.a,  parm.b,  parm.c))
    
    
    # Sp2 Oxygen
    parm = parms['O']['Sp2']
    print('\n\n{:^24}'.format('Sp2 Oxygen'))
    print('{:^8} {:^8} {:^8}'.format('a', 'b', 'c'))
    print('{:^8} {:^8} {:^8}'.format(parm.a,  parm.b,  parm.c))
    
    
    # Hydrogen
    parm = parms['H']['Terminal']
    print('\n\n{:^24}'.format('Hydrogen'))
    print('{:^8} {:^8} {:^8}'.format('a', 'b', 'c'))
    print('{:^8} {:^8} {:^8}'.format(parm.a,  parm.b,  parm.c))
    
    
    # Sulfur
    parm = parms['S']['Sp3']
    print('\n\n{:^24}'.format('Sp3 Sulfur'))
    print('{:^8} {:^8} {:^8}'.format('a', 'b', 'c'))
    print('{:^8} {:^8} {:^8}'.format(parm.a,  parm.b,  parm.c))

