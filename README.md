# LUNAR (version 05September2025)
LUNAR stands for LAMMPS Utility (for) Network Analysis (and) Reactivity and is a stand alone Python (3.7+) toolkit to supplement LAMMPS. LUNAR is focused on pre-processing and post-processing inputs and outputs of LAMMPS with emphasis of using LAMMPS for producing structure-property relationships for ICME process modeling of polymers. However many of the codes within LUNAR can be used outside of process modeling type of Molecular Dynamics simulations.

LUNAR was written by Josh Kemppainen during his PhD at Michigan Technological Univeristy with Dr. Gregory M. Odegard being his research advisor. LUNAR is maintained by Josh Kemppainen and Dr. Jacob R. Gissinger. LUNAR is not planned on being packaged with PyPi, since all the source code is provided and allows the user to modify the source code easily. If users of LUNAR add considerably substantial modules please notify Josh Kemppainen or Dr. Jacob R. Gissinger to review the additions, such that they maybe included in the LUNAR distribution and maintained for future releases.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

## The LUNAR distribution includes the following files and directories:
- **README**                     this file
- **LICENSE**                    the GNU General Public License (GPL)
- **DEV_Tools**                  undocumented codes that may be of use to others
- **doc**                        documentation
- **EXAMPLES**                   simple test problems
- **frc_files**                  included interatomic potential files
- **src**                        source files
- **all2lmp.py**                 LUNAR file to run all2lmp.py
- **atom_removal.py**            LUNAR file to run atom_removal.py
- **atom_typing.py**             LUNAR file to run atom_typing.py
- **auto_cluster_analysis.py**   LUNAR file to run auto_cluster_analysis.py
- **auto_morse_bond_update.py**  LUNAR file to run all2lmp.py
- **bond_react_merge.py**        LUNAR file to run bond_react_merge.py
- **bond_react_merge_prep.py**   LUNAR file to run bond_react_merge_prep.py
- **cell_builder.py**            LUNAR file to run cell_builder.py
- **cluster_analysis.py**        LUNAR file to run cluster_analysis.py
- **convert2graphite.py**        LUNAR file to run convert2graphite.py
- **free_volume.py**             LUNAR file to run free_volume.py
- **LUNAR.py**                   LUNAR file to run LUNAR.py
- **lmp2SYBYLmol2.py**           LUNAR file to run lmp2SYBYLmol2.py
- **log_analysis.py**            LUNAR file to run log_analysis.py
- **sheet_builder.py**           LUNAR file to run sheet_builder.py

## Installation
- Install python (3.7+) via a number of different methods such as:
  1. https://www.python.org/
  2. https://www.anaconda.com/
  3. sudo apt-get install python3.11 (Linux Ubuntu example)

- LUNAR will run on the standard python library where some of the code/features requires the following third party modules to be installed (preferable via pip):
  - matplotlib
  - numpy
  - numba
  - tqdm
  - rdkit
  - natsort
  - scipy
  - pwlf
- Please resort to the documention of your Python installion method on the syntax for using pip before reaching out to others to help install the third party libraries listed above. In general you will use pip as:
  - https://www.python.org/
    - python -m pip install PACKAGE
  - https://www.anaconda.com/
    - pip install PACKAGE
  - sudo apt-get install python3.11 (Linux Ubuntu example)
    - pip install PACKAGE

- Note that Python v3.11+ will run about 1.5x as quick vs v3.10, thus it is recommended to install v.3.11+ even though the code assumes v3.7+ features.

## Terms of use
- You agree to acknowledge the use of the LUNAR in any reports or publications of results obtained with the LUNAR by citing "LUNAR: Automated Input Generation and Analysis for Reactive LAMMPS Simulations" at https://doi.org/10.1021/acs.jcim.4c00730.
- LUNAR is not responsible for the validity contents of any linked website and does not accept any liability for actions or consequences applying LUNAR codes.

## Additional resources
- Tutorials and workshops related to LUNAR and LAMMPS will be posted at: https://www.youtube.com/@CMMRLab-LUNAR.
- Files generated during posted tutorials on YouTube will be posted at: https://github.com/CMMRLab/LUNAR_tutorials/tree/main.
- Feel free to add comments on the YouTube channel of questions.
