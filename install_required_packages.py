# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.1
April 23rd, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

    *****************************************
    * Requirements:                         *
    *   python 3.7+                         *
    *                                       *
    * Run mode as:                          *
    *   python install_python_packages.py   *
    *                or                     *
    *   python3 install_python_packages.py  *
    *                                       *
    *   Depending on how python is aliased  *
    *                                       *
    * NOTE that if using anaconda Spyder,   *
    * this code can NOT be run from the IDE *
    * and must be run through Anaconda      *
    * prompt                                *
    *****************************************
"""
# Import Necessary Libraries
import subprocess
import sys


# Required packages
required_packages = ['matplotlib', 'numpy', 'numba', 'tqdm', 'rdkit', 'natsort', 'scipy', 'pwlf']
print(f'Will attempt to install the following packages: {", ". join(required_packages)}\n\n')


# Function to install packages
def install(package):
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', package])
    
    
# Install packages
installed = []; failed = [];
for package in required_packages:
    try: 
        install(package)
        installed.append(installed)
    except: 
        failed.append(package)
        print(f'ERROR {package} failed to install, resort to manual install methods.')
        
        
# Print results of install
print(f'\n\nInstalled {len(installed)} of {len(required_packages)} packages.')
if failed:
    for package in failed:
        print(f'{package} failed to be installed, resort to manual install methods.')
    