# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
September 6th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Function to find setting to update and return new line with updated setting
def parse_and_modify(line, setting, stringflag=False, splitchar='='):
    newline = '';  splitline = line.split(splitchar);
    newline += splitline[0]
    newline += '{} '.format(splitchar)
    if stringflag: newline += "'{}'".format(setting)
    else: newline += str(setting)
    commentflag = False; commaflag = False;
    ending = ''.join([i for n, i in enumerate(splitline) if n > 0])
    for i in ending:
        if i in [',', '{', '}'] and splitchar == ':':
            commaflag = True
        if i == '#':
            commentflag = True
            if not commaflag: newline += ' '
        if commentflag or commaflag:
            newline += i
    if not commentflag:
        newline += '\n'
    return newline

# Function to read pyscript
def read(pyscript):
    lines = []
    with open(pyscript, 'r') as f:
        for wholeline in f:
            lines.append(wholeline)
    return lines