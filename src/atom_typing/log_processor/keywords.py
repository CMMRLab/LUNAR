# -*- coding: utf-8 -*-
"""
Revision 1.0
February 20, 2025
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""

# Set up a map of shortcut strings to keyword shortcuts.
# shortcuts={'shortcut-string': 'keyword-string'}
shortcuts = {'file': 'filename',
             
             '%Ring-Count:3': "rings[3]['%Ring']",
             '%Ring-Count:4': "rings[4]['%Ring']",
             '%Ring-Count:5': "rings[5]['%Ring']",
             '%Ring-Count:6': "rings[6]['%Ring']",
             '%Ring-Count:7': "rings[7]['%Ring']",
             
             '%Ring-Mass:C3': "rings[6]['C']['%Mass']",
             '%Ring-Mass:C4': "rings[4]['C']['%Mass']",
             '%Ring-Mass:C5': "rings[5]['C']['%Mass']",
             '%Ring-Mass:C6': "rings[6]['C']['%Mass']",
             '%Ring-Mass:C7': "rings[7]['C']['%Mass']",

             '%Ring-Count:C3': "rings[6]['C']['%natoms']",
             '%Ring-Count:C4': "rings[4]['C']['%natoms']",
             '%Ring-Count:C5': "rings[5]['C']['%natoms']",
             '%Ring-Count:C6': "rings[6]['C']['%natoms']",
             '%Ring-Count:C7': "rings[7]['C']['%natoms']",

             '%Hybrid-Mass:all-H': "hybridizations['all-H']['%Mass']",
             '%Hybrid-Count:all-H': "hybridizations['all-H']['%natoms']",
             
             '%Hybrid-Mass:Sp1-C': "hybridizations['Sp1-C']['%Mass']",
             '%Hybrid-Mass:Sp2-C': "hybridizations['Sp2-C']['%Mass']",
             '%Hybrid-Mass:Sp3-C': "hybridizations['Sp3-C']['%Mass']",
             '%Hybrid-Mass:all-C': "hybridizations['all-C']['%Mass']",
             '%Hybrid-Count:Sp1-C': "hybridizations['Sp1-C']['%natoms']",
             '%Hybrid-Count:Sp2-C': "hybridizations['Sp2-C']['%natoms']",
             '%Hybrid-Count:Sp3-C': "hybridizations['Sp3-C']['%natoms']",
             '%Hybrid-Count:all-C': "hybridizations['all-C']['%natoms']",
             
             '%Hybrid-Mass:Sp1-O': "hybridizations['Sp1-O']['%Mass']",
             '%Hybrid-Mass:Sp2-O': "hybridizations['Sp2-O']['%Mass']",
             '%Hybrid-Mass:Sp3-O': "hybridizations['Sp3-O']['%Mass']",
             '%Hybrid-Mass:all-O': "hybridizations['all-O']['%Mass']",
             '%Hybrid-Count:Sp1-O': "hybridizations['Sp1-O']['%natoms']",
             '%Hybrid-Count:Sp2-O': "hybridizations['Sp2-O']['%natoms']",
             '%Hybrid-Count:Sp3-O': "hybridizations['Sp3-O']['%natoms']",
             '%Hybrid-Count:all-O': "hybridizations['all-O']['%natoms']",
             
             '%Hybrid-Mass:Sp1-N': "hybridizations['Sp1-N']['%Mass']",
             '%Hybrid-Mass:Sp2-N': "hybridizations['Sp2-N']['%Mass']",
             '%Hybrid-Mass:Sp3-N': "hybridizations['Sp3-N']['%Mass']",
             '%Hybrid-Mass:all-N': "hybridizations['all-N']['%Mass']",
             '%Hybrid-Count:Sp1-N': "hybridizations['Sp1-N']['%natoms']",
             '%Hybrid-Count:Sp2-N': "hybridizations['Sp2-N']['%natoms']",
             '%Hybrid-Count:Sp3-N': "hybridizations['Sp3-N']['%natoms']",
             '%Hybrid-Count:all-N': "hybridizations['all-N']['%natoms']",
             }